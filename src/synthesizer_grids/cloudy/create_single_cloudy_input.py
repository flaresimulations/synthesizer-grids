"""
Create a single cloudy script. This is used for testing.
"""

import numpy as np
import argparse
from pathlib import Path
import yaml
import h5py

# synthesiser modules
from synthesizer.abundances import (
    Abundances,
    depletion_models,
    solar_abundance_patterns,
)

from synthesizer.grid import Grid
from synthesizer.photoionisation import cloudy17, cloudy23


# local modules
from utils import get_grid_properties, apollo_submission_script


def load_grid_params(param_file="c17.03-sps", param_dir="params"):
    """
    Read parameters from a yaml parameter file

    Arguments:
        param_file (str)
            filename of the parameter file
        param_dir (str)
            directory containing the parameter file

    Returns:
        fixed_params (dict)
            dictionary of parameters that are fixed
        grid_params (dict)
            dictionary of parameters that vary on the grid
    """

    # open paramter file
    with open(f"{param_dir}/{param_file}.yaml", "r") as stream:
        try:
            params = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    grid_params = {}
    fixed_params = {}

    # loop over parameters
    for k, v in params.items():

        # if parameter is a list store it in the grid_parameters dictionary
        # and convert to a numpy array
        if isinstance(v, list):
            grid_params[k] = np.array(list(map(float, v)))

        # otherwise store it in fixed_params dictionary
        else:
            fixed_params[k] = v

    return fixed_params, grid_params


if __name__ == "__main__":

    machine = 'apollo'
    incident_grid = 'bc03-2016-BaSeL_chabrier-0.1,100'
    grid_dir = '/Users/sw376/Dropbox/Research/data/synthesizer/grids'
    cloudy_dir = '/Users/sw376/Dropbox/Research/data/synthesizer/cloudy'
    cloudy_params = 'c23.01-sps-c17.03'
    cloudy_params = 'c17.03-sps'

    verbose = True

    # load the cloudy parameters you are going to run
    fixed_params, grid_params = load_grid_params(cloudy_params)

    # set cloudy version
    if fixed_params['cloudy_version'] == 'c17.03':
        cloudy = cloudy17
    if fixed_params['cloudy_version'] == 'c23.01':
        cloudy = cloudy23
    print(cloudy)




    # open the parent incident grid
    incident_grid = Grid(
        incident_grid,
        grid_dir=f"{grid_dir}",
        read_lines=False,
    )

    # get name of new grid (concatenation of incident_grid and cloudy parameter file)
    new_grid_name = f"{incident_grid}_cloudy-{cloudy_params}"

    # set a list of the axes
    axes = list(incident_grid.axes) + list(grid_params.keys())
    if verbose:
        print("axes:", axes)

    # add the incident grid parameters to grid_params
    for axis in incident_grid.axes:
        grid_params[axis] = getattr(incident_grid, axis)

    if verbose:
        # print fixed parameters
        for k, v in fixed_params.items():
            print(k, v)

        # print grid parameters, including incident parameters
        for k, v in grid_params.items():
            print(k, v)

    # if the U model is the reference model (i.e. not fixed) save the grid point for the reference values
    if fixed_params["ionisation_parameter_model"] == "ref":
        # get the indices of the reference grid point (this is used by the reference model)
        incident_ref_grid_point = incident_grid.get_grid_point(
            [fixed_params["reference_" + k] for k in incident_grid.axes]
        )

        # add these to the parameter file
        for k, i in zip(incident_grid.axes, incident_ref_grid_point):
            fixed_params["reference_" + k + "_index"] = i

    # combine all parameters
    params = fixed_params | grid_params

    # get properties of the grid
    (
        n_axes,
        shape,
        n_models,
        mesh,
        model_list,
        index_list,
    ) = get_grid_properties(axes, grid_params, verbose=True)

    
    # loop over all models
    for i, (grid_params_tuple, grid_index_tuple) in enumerate(
        zip(model_list[:1], index_list[:1])
    ):
        # get a dictionary of all parameters
        grid_params_ = dict(zip(axes, grid_params_tuple))

        # get a dictionary of the parameter grid point
        grid_index_ = dict(zip(axes, grid_index_tuple))

        # get a dictionary of just the incident parameters
        incident_params_ = {k: grid_params_[k] for k in incident_grid.axes}

        # get a dictionary of the incident parameter grid point
        incident_index_ = {k: grid_index_[k] for k in incident_grid.axes}

        # get a tuple of the incident grid point
        incident_grid_point = tuple(grid_index_[k] for k in incident_grid.axes)

        # join the fixed and current iteration of the grid parameters
        params_ = fixed_params | grid_params_

        # set cloudy metallicity parameter to the stellar metallicity
        if "metallicity" in grid_params_.keys():
            params_["metallicity"] = grid_params_["metallicity"]
        elif "log10metallicity" in grid_params_.keys():
            params_["metallicity"] = 10 ** grid_params_["log10metallicity"]

        # create abundance object
        abundances = Abundances(
            metallicity=float(params_["metallicity"]),
            solar=params_["solar_abundance"],
            alpha=params_["alpha"],
            abundances=params_["abundance_scalings"],
            depletion_model=params_["depletion_model"],
            depletion_scale=params_["depletion_scale"],
        )

        # if reference U model is used
        if params_["ionisation_parameter_model"] == "ref":
            # Calculate the difference between the reference 
            # log10_specific_ionising_luminosity for HI and the current grid
            # point
            delta_log10_specific_ionising_luminosity = (
                incident_grid.log10_specific_ionising_lum["HI"][
                    incident_grid_point
                ]
                - incident_grid.log10_specific_ionising_lum["HI"][
                    incident_ref_grid_point
                ]
            )

            # for spherical geometry the effective log10U is this
            if params_["geometry"] == "spherical":
                log10U = (
                    np.log10(params_["reference_ionisation_parameter"])
                    + (1 / 3) * delta_log10_specific_ionising_luminosity
                )

            # For plane-parallel geometry the effective just scales with
            # log10_specific_ionising_luminosity
            elif params_["geometry"] == "planeparallel":
                log10U = (
                    np.log10(params_["reference_ionisation_parameter"])
                    + delta_log10_specific_ionising_luminosity
                )

            else:
                print(
                    f"ERROR: do not understand geometry choice: {params_['geometry']}"
                )

        # if fixed U model is used
        elif params_["ionisation_parameter_model"] == "fixed":
            log10U = np.log10(params_["ionisation_parameter"])

        else:
            print(
                f"ERROR: do not understand U model choice: {params_['ionisation_parameter_model']}"
            )

        # set log10U to provide cloudy
        params_["ionisation_parameter"] = 10 ** float(log10U)

        # get wavelength
        lam = incident_grid.lam  # AA

        # get luminosity
        lnu = incident_grid.spectra["incident"][incident_grid_point]

        # this returns the relevant shape commands, in this case for a tabulated SED
        shape_commands = cloudy.ShapeCommands.table_sed(
            str(i + 1), lam, lnu, output_dir='test/'
        )

        # create cloudy input file
        cinput = cloudy.create_cloudy_input(
            str(i + 1),
            shape_commands,
            abundances,
            output_dir='test/',
            **params_,
        )

        
        print(''.join(cinput))
        
  