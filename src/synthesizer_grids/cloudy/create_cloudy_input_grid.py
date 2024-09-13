"""
Create a grid of cloudy scripts based on combining an incident grid with a set
of cloudy parameters. Also creates a machine specific script.
"""

from pathlib import Path

import numpy as np
import yaml
from synthesizer.abundances import (
    Abundances,
)
from synthesizer.grid import Grid
from synthesizer.photoionisation import cloudy17, cloudy23

# local modules
from utils import apollo_submission_script, get_grid_properties

from synthesizer_grids.parser import Parser


def load_grid_params(param_file="c23.01-sps", param_dir="params"):
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

    # Loop over parameters
    for k, v in params.items():
        # If parameter is a list store it in the grid_parameters dictionary
        # and convert to a numpy array
        if isinstance(v, list):
            grid_params[k] = np.array(list(map(float, v)))

        # Otherwise store it in fixed_params dictionary
        else:
            fixed_params[k] = v

    return fixed_params, grid_params


if __name__ == "__main__":
    parser = Parser(description="Run a grid of incident cloudy models")

    # Machine (for submission script generation)
    parser.add_argument("-machine", type=str, required=True)

    # Path to directory where cloudy runs are
    parser.add_argument("-cloudy_dir", type=str, required=True)

    # The name of the incident grid
    parser.add_argument("-incident_grid", type=str, required=True)

    # The cloudy reference parameter set
    parser.add_argument(
        "-cloudy_params", type=str, required=True, default="c23.01-sps"
    )

    # A second cloudy parameter set which supersedes the above
    parser.add_argument(
        "-cloudy_params_addition",
        type=str,
        required=False,
    )

    # Path to cloudy directory (not the executable; this is assumed to
    # {cloudy}/{cloudy_version}/source/cloudy.exe)
    parser.add_argument("-cloudy_path", type=str, required=False)

    # Verbosity flag
    parser.add_argument("-verbose", type=bool, required=False, default=True)

    # Parse arguments
    args = parser.parse_args()

    verbose = args.verbose

    # Load the cloudy parameters you are going to run
    fixed_params, grid_params = load_grid_params(args.cloudy_params)

    # If an additional parameter set is provided supersede the default
    # parameters with these.
    if args.cloudy_params_addition:
        additional_fixed_params, additional_grid_params = load_grid_params(
            args.cloudy_params_addition
        )
        fixed_params = fixed_params | additional_fixed_params
        grid_params = grid_params | additional_grid_params

    # Set cloudy version
    if fixed_params["cloudy_version"] == "c17.03":
        cloudy = cloudy17
    if fixed_params["cloudy_version"] == "c23.01":
        cloudy = cloudy23
    print(cloudy)

    # Open the parent incident grid
    incident_grid = Grid(
        args.incident_grid + ".hdf5",
        grid_dir=f"{args.grid_dir}",
        read_lines=False,
    )

    # Get name of new grid (concatenation of incident_grid and cloudy
    # parameter file)
    new_grid_name = f"{args.incident_grid}_cloudy-{args.cloudy_params}"

    # If an additional parameter set append this to the new grid name
    if args.cloudy_params_addition:
        # ignore the directory part
        cloudy_params_addition_name = args.cloudy_params_addition.split("/")[
            -1
        ]
        new_grid_name += "-" + cloudy_params_addition_name

    # Define output directories
    output_dir = f"{args.cloudy_dir}/{new_grid_name}"

    # Make output directories
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # For submission system output files
    Path(f"{output_dir}/output").mkdir(parents=True, exist_ok=True)

    # Set a list of the axes
    axes = list(incident_grid.axes) + list(grid_params.keys())
    if verbose:
        print("axes:", axes)

    # Add the incident grid parameters to grid_params
    for axis in incident_grid.axes:
        grid_params[axis] = getattr(incident_grid, axis)

    if verbose:
        # print fixed parameters
        for k, v in fixed_params.items():
            print(k, v)

        # print grid parameters, including incident parameters
        for k, v in grid_params.items():
            print(k, v)

    # If the ionisation_parameter_model is the reference model (i.e. not fixed)
    # save the grid point for the reference values
    if fixed_params["ionisation_parameter_model"] == "ref":
        # get the indices of the reference grid point (this is used by the
        # reference model)
        print(
            "grid point to get:",
            type(fixed_params["reference_ages"]),
            type(fixed_params["reference_metallicities"]),
        )
        incident_ref_grid_point = incident_grid.get_grid_point(
            [fixed_params["reference_" + k] for k in incident_grid.axes]
        )

        # add these to the parameter file
        for k, i in zip(incident_grid.axes, incident_ref_grid_point):
            fixed_params["reference_" + k + "_index"] = i

    # Combine all parameters
    params = fixed_params | grid_params

    # Add the parameter file as a parameter
    params["parameter_file"] = args.cloudy_params

    # Save all parameters
    yaml.dump(params, open(f"{output_dir}/params.yaml", "w"))

    # Get properties of the grid
    (
        n_axes,
        shape,
        n_models,
        mesh,
        model_list,
        index_list,
    ) = get_grid_properties(axes, grid_params, verbose=True)

    # Loop over all models
    for i, (grid_params_tuple, grid_index_tuple) in enumerate(
        zip(model_list, index_list)
    ):
        # Get a dictionary of all parameters
        grid_params_ = dict(zip(axes, grid_params_tuple))

        # Get a dictionary of the parameter grid point
        grid_index_ = dict(zip(axes, grid_index_tuple))

        # Get a dictionary of just the incident parameters
        incident_params_ = {k: grid_params_[k] for k in incident_grid.axes}

        # Get a dictionary of the incident parameter grid point
        incident_index_ = {k: grid_index_[k] for k in incident_grid.axes}

        # Get a tuple of the incident grid point
        incident_grid_point = tuple(grid_index_[k] for k in incident_grid.axes)

        # Join the fixed and current iteration of the grid parameters
        params_ = fixed_params | grid_params_

        # Set cloudy metallicities parameter to the stellar metallicities
        if "metallicities" in grid_params_.keys():
            params_["metallicities"] = grid_params_["metallicities"]
        elif "log10metallicities" in grid_params_.keys():
            params_["metallicities"] = 10 ** grid_params_["log10metallicities"]

        # Create abundance object
        abundances = Abundances(
            metallicity=float(params_["metallicities"]),
            reference=params_["reference_abundance"],
            alpha=params_["alpha"],
            abundances=params_["abundance_scalings"],
            depletion_model=params_["depletion_model"],
            depletion_scale=params_["depletion_scale"],
        )

        # If reference U model is used
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

            # For spherical geometry the effective log10U is this
            if params_["geometry"] == "spherical":
                log10U = (
                    np.log10(params_["reference_ionisation_parameter"])
                    + (1 / 3) * delta_log10_specific_ionising_luminosity
                )

            # For spherical geometry the effective log10U is this.
            # The spherical-U model uses U as a direct input to cloudy
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

        # If fixed U model is used:
        elif params_["ionisation_parameter_model"] == "fixed":
            log10U = np.log10(params_["ionisation_parameter"])

        else:
            print(
                f"""ERROR: do not understand U model choice:
                {params_['ionisation_parameter_model']}"""
            )

        # Set log10U to provide cloudy
        params_["ionisation_parameter"] = 10 ** float(log10U)

        # Get wavelength
        lam = incident_grid.lam  # AA

        # Get luminosity
        lnu = incident_grid.spectra["incident"][incident_grid_point]

        # This returns the relevant shape commands, in this case for a
        # tabulated SED
        shape_commands = cloudy.ShapeCommands.table_sed(
            str(i + 1), lam, lnu, output_dir=output_dir
        )

        # Create cloudy input file
        cloudy.create_cloudy_input(
            str(i + 1),
            shape_commands,
            abundances,
            output_dir=output_dir,
            **params_,
        )

    # Create submission script
    if args.machine == "apollo":
        apollo_submission_script(
            n_models, output_dir, args.cloudy_path, params_["cloudy_version"]
        )
