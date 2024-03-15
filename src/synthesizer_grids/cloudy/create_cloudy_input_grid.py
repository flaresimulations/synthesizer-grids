"""
Create a grid of cloudy scripts based on combining an incident grid with a set
of cloudy parameters. Also creates a machine specific script.
"""

import numpy as np
import argparse
from pathlib import Path
import yaml
import h5py

from synthesizer.abundances import (
    Abundances,
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
    parser = argparse.ArgumentParser(
        description="Run a grid of incident cloudy models"
    )

    # machine (for submission script generation)
    parser.add_argument("-machine", type=str, required=True)

    # path to grid directory (i.e. where incident and new grids are stored)
    parser.add_argument("-grid_dir", type=str, required=True)

    # path to directory where cloudy runs are
    parser.add_argument("-cloudy_dir", type=str, required=True)

    # the name of the incident grid
    parser.add_argument("-incident_grid", type=str, required=True)

    # the cloudy parameters, including any grid axes
    parser.add_argument(
        "-cloudy_params",
        type=str,
        required=False,
        default="c17.03-sps"
    )

    # path to cloudy directory (not executable; this is assumed to
    # {cloudy}/{cloudy_version}/source/cloudy.exe)
    parser.add_argument("-cloudy_path", type=str, required=False)

    # verbosity flag
    parser.add_argument("-verbose", type=bool, required=False, default=True)

    # parse arguments
    args = parser.parse_args()

    verbose = args.verbose

    # load the cloudy parameters you are going to run
    fixed_params, grid_params = load_grid_params(args.cloudy_params)

    # set cloudy version
    if fixed_params['cloudy_version'] == 'c17.03':
        cloudy = cloudy17
    if fixed_params['cloudy_version'] == 'c23.01':
        cloudy = cloudy23
    print(cloudy)

    # open the parent incident grid
    incident_grid = Grid(
        args.incident_grid,
        grid_dir=f"{args.grid_dir}",
        read_lines=False,
    )

    # get name of new grid (concatenation of incident_grid and cloudy
    # parameter file)
    new_grid_name = f"{args.incident_grid}_cloudy-{args.cloudy_params}"

    # define output directories
    output_dir = f"{args.cloudy_dir}/{new_grid_name}"

    # make output directories
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # for submission system output files
    Path(f"{output_dir}/output").mkdir(parents=True, exist_ok=True)

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

    # if the ionisation_parameter_model is the reference model (i.e. not fixed)
    # save the grid point for the reference values
    if fixed_params["ionisation_parameter_model"] == "ref":
        # get the indices of the reference grid point (this is used by the
        # reference model)
        incident_ref_grid_point = incident_grid.get_grid_point(
            [fixed_params["reference_" + k] for k in incident_grid.axes]
        )

        # add these to the parameter file
        for k, i in zip(incident_grid.axes, incident_ref_grid_point):
            fixed_params["reference_" + k + "_index"] = i

    # combine all parameters
    params = fixed_params | grid_params

    # add the parameter file as a parameter
    params['parameter_file'] = args.cloudy_params

    # save all parameters
    yaml.dump(params, open(f"{output_dir}/params.yaml", "w"))

    # get properties of the grid
    (
        n_axes,
        shape,
        n_models,
        mesh,
        model_list,
        index_list,
    ) = get_grid_properties(axes, grid_params, verbose=True)

    # create new synthesizer grid to contain the new grid
    # open the new grid
    with h5py.File(f"{args.grid_dir}/{new_grid_name}.hdf5", "w") as hf:
        # open the original incident model grid
        with h5py.File(
            f"{args.grid_dir}/{args.incident_grid}.hdf5", "r"
        ) as hf_incident:

            # print out datasets
            if verbose:
                hf_incident.visit(print)

            # copy top-level attributes
            for k, v in hf_incident.attrs.items():

                # If v is None then convert to string None for saving in the 
                # HDF5 file.
                if v is None:
                    v = 'None'

                hf.attrs[k] = v
                if verbose:
                    print(k, v)

            # add attribute with the original incident grid axes
            hf.attrs["incident_axes"] = hf_incident.attrs["axes"]

            # we want to copy over log10_specific_ionising_luminosity from the
            # incident grid to allow us to normalise the cloudy outputs.
            # However, the axes of the incident grid may be different from the
            # cloudy grid due to additional parameters, in which we need to
            # extend the axes of log10_specific_ionising_luminosity.

            # if there are no additional axes simply copy over the incident
            # log10_specific_ionising_luminosity
            if len(axes) == len(hf.attrs["incident_axes"]):
                hf_incident.copy("log10_specific_ionising_luminosity", hf)

            # else we need to expand the axis
            else:
                # this is amount by which we need to expand
                expansion = int(
                    np.product(shape)
                    / np.product(
                        hf_incident[
                            "log10_specific_ionising_luminosity/HI"
                        ].shape
                    )
                )

                # loop over ions
                for ion in hf_incident[
                    "log10_specific_ionising_luminosity"
                ].keys():
                    # get the incident log10_specific_ionising_luminosity array
                    log10_specific_ionising_luminosity_incident = hf_incident[
                        f"log10_specific_ionising_luminosity/{ion}"
                    ][()]

                    # create new array with repeated elements
                    log10_specific_ionising_luminosity = np.repeat(
                        log10_specific_ionising_luminosity_incident,
                        expansion,
                        axis=-1,
                    )

                    # reshape array to match new shape and save
                    hf[
                        f"log10_specific_ionising_luminosity/{ion}"
                    ] = np.reshape(log10_specific_ionising_luminosity, shape)

        # add attribute with full grid axes
        hf.attrs["axes"] = axes

        # add the bin centres for the grid bins
        for axis in axes:
            hf[f"axes/{axis}"] = grid_params[axis]

        # add other parameters as attributes
        for k, v in params.items():

            # If v is None then convert to string None for saving in the 
            # HDF5 file.
            if v is None:
                v = 'None'

            # if the parameter is a dictionary (e.g. as used for abundances)
            if isinstance(v, dict):
                for k2, v2 in v.items():
                    hf.attrs[k+'_'+k2] = v2
            else:
                hf.attrs[k] = v

        if verbose:
            print("-" * 50)
            print("---- attributes")
            for k, v in hf.attrs.items():
                print(k, v)
            print("---- groups and datasets")
            hf.visit(print)

    # loop over all models
    for i, (grid_params_tuple, grid_index_tuple) in enumerate(
        zip(model_list, index_list)
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
            reference=params_["reference_abundance"],
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

            # for spherical geometry the effective log10U is this.
            # the spherical-U model uses U as a direct input to cloudy 
            # instead of calculating Q.
            elif params_["geometry"] == "spherical-U":
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
                    f"""ERROR: do not understand geometry choice:
                    {params_['geometry']}"""
                )

        # if fixed U model is used
        elif params_["ionisation_parameter_model"] == "fixed":
            log10U = np.log10(params_["ionisation_parameter"])

        else:
            print(
                f"""ERROR: do not understand U model choice:
                {params_['ionisation_parameter_model']}"""
            )

        # set log10U to provide cloudy
        params_["ionisation_parameter"] = 10 ** float(log10U)

        # get wavelength
        lam = incident_grid.lam  # AA

        # get luminosity
        lnu = incident_grid.spectra["incident"][incident_grid_point]

        # this returns the relevant shape commands, in this case for a
        # tabulated SED
        shape_commands = cloudy.ShapeCommands.table_sed(
            str(i + 1), lam, lnu, output_dir=output_dir
        )

        # create cloudy input file
        cloudy.create_cloudy_input(
            str(i + 1),
            shape_commands,
            abundances,
            output_dir=output_dir,
            **params_,
        )

    # create submission script
    if args.machine == "apollo":
        apollo_submission_script(
            n_models, output_dir, args.cloudy_path, params_["cloudy_version"]
        )
