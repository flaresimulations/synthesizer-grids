"""A module for reading in cloudy outputs and creating a synthesizer grid.

This module reads in cloudy outputs and creates a synthesizer grid. The
synthesizer grid is a HDF5 file that contains the spectra and line luminosities
from cloudy outputs.

Example usage:
    python create_synthesizer_grid.py --cloudy_dir=cloudy_dir
    --incident_grid=incident_grid --cloudy_params=cloudy_params
"""

import os
import shutil

import h5py
import numpy as np
import yaml
from synthesizer.grid import Grid
from synthesizer.photoionisation import cloudy17, cloudy23
from synthesizer.sed import Sed
from synthesizer.units import has_units
from unyt import Angstrom, Hz, cm, dimensionless, erg, eV, s, yr
from utils import get_grid_props_cloudy

from synthesizer_grids.grid_io import GridFile
from synthesizer_grids.parser import Parser

# Create dictionary of commonly used units for grid axes
axes_units = {
    "reference_ionisation_parameter": dimensionless,
    "ionisation_parameter": dimensionless,
    "ages": yr,
    "metallicities": dimensionless,
    "hydrogen_density": cm ** (-3),
}


def create_empty_grid(
    grid_dir,
    incident_grid_name,
    new_grid_name,
    params,
    grid_params,
):
    """
    Create an empty Synthesizer grid file based on an existing incident grid.

    Args:
        grid_dir (str):
            Directory where the grid file will be saved.
        incident_grid_name (str):
            Name of the incident grid file we ran through cloudy.
        new_grid_name (str):
            Name of the new grid file.
        params (dict):
            Dictionary of parameters to be saved in the grid file.
        grid_params (dict):
            Dictionary of grid parameters.
    """
    # Open the parent incident grid
    incident_grid = Grid(
        incident_grid_name + ".hdf5",
        grid_dir=grid_dir,
        read_lines=False,
    )

    # Make the new grid
    out_grid = GridFile(f"{grid_dir}/{new_grid_name}.hdf5")

    # Add attribute with the original incident grid axes
    out_grid.write_attribute(
        group="/",
        attr_key="incident_axes",
        data=incident_grid.axes,
    )

    # Set a list of the axes including the incident and new grid axes
    axes = list(incident_grid.axes) + list(grid_params.keys())

    # Add the incident grid parameters to grid_params
    for axis in incident_grid.axes:
        grid_params[axis] = getattr(incident_grid, axis)

    # Get properties of the grid
    (
        n_axes,
        shape,
        n_models,
        mesh,
        model_list,
        index_list,
    ) = get_grid_props_cloudy(axes, grid_params, verbose=True)

    # We want to copy over log10_specific_ionising_luminosity from the
    # incident grid to allow us to normalise the cloudy outputs.
    # However, the axes of the incident grid may be different from the
    # cloudy grid due to additional parameters, in which case we need to
    # extend the axes of log10_specific_ionising_luminosity.

    # Compute the amount we may have to expand the axes (if we don't need to
    # expand this is harmless)
    expansion = int(
        np.product(shape)
        / np.product(
            incident_grid.log10_specific_ionising_luminosity["HI"].shape
        )
    )

    # Loop over ions
    for ion in incident_grid.log10_specific_ionising_luminosity.keys():
        # If there are no additional axes simply copy over the incident
        # log10_specific_ionising_luminosity .
        if len(axes) == len(incident_grid.axes):
            out_grid.write_dataset(
                key=f"log10_specific_ionising_luminosity/{ion}",
                data=incident_grid.log10_specific_ionising_luminosity[ion]
                * dimensionless,
                description="The specific ionising photon luminosity for"
                f" {ion}",
                log_on_read=False,
            )

        # Otherwise, we need to expand the axis to account for these additional
        # parameters.
        else:
            # Create new array with repeated elements
            log10_specific_ionising_luminosity = np.repeat(
                incident_grid.log10_specific_ionising_luminosity[ion]
                * dimensionless,
                expansion,
                axis=-1,
            )
            out_grid.write_dataset(
                f"log10_specific_ionising_luminosity/{ion}",
                np.reshape(log10_specific_ionising_luminosity, shape)
                * dimensionless,
                description="The specific ionising photon luminosity"
                f" for {ion}",
                log_on_read=False,
            )

    # Add attribute with the full grid axes
    out_grid.write_attribute(group="/", attr_key="axes", data=axes)

    # Now write out the grid axes, we first do the incident grid axes so we
    # can extract their metadata and then any extras
    for axis, log_on_read in zip(
        incident_grid.axes, incident_grid._logged_axes
    ):
        # Write out this axis (we can get everything we need from the incident)
        out_grid.write_dataset(
            key=f"axes/{axis}",
            data=getattr(incident_grid, axis),
            description=f"Grid axes {axis} (taken from incident grid)",
            log_on_read=log_on_read,
        )

    # Now write out the new grid axes if there are any
    for axis in axes:
        # Skip the one we did above
        if axis in incident_grid.axes:
            continue

        # If we have units in the grid_params dictionary already then use
        # these, otherwise use the axes_units dictionary, if we don't have
        # units for the axis then raise an error.
        if has_units(grid_params[axis]):
            out_grid.write_dataset(
                key=f"axes/{axis}",
                data=grid_params[axis],
                description=f"grid axes {axis} (introduced during cloudy run)",
                log_on_read=False,
            )

        # Ok, maybe its just a common axis that we have units for
        elif axis in axes_units:
            # Multiply values with the corresponding unit
            values_with_units = grid_params[axis] * axes_units[axis]

            out_grid.write_dataset(
                key=f"axes/{axis}",
                data=values_with_units,
                description=f"grid axes {axis} (introduced during cloudy run)",
                log_on_read=False,
            )

        else:
            raise ValueError(
                f"Unit for axis '{axis}' needs to be added to "
                f"create_synthesizer_grid.py"
            )

    # Add other parameters as attributes
    for k, v in params.items():
        # If v is None then convert to string None for saving in the
        # HDF5 file.
        if v is None:
            v = "None"

        # If the parameter is a dictionary (e.g. as used for abundances)
        if isinstance(v, dict):
            for k2, v2 in v.items():
                out_grid.write_attribute(group="k", attr_key=k2, data=v2)
        else:
            out_grid.write_attribute(group="/", attr_key=k, data=v)


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

    # open parameter file
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


def get_grid_properties_hf(hf, verbose=True):
    """
    A wrapper over get_grid_properties to get the grid properties for a HDF5
    grid.
    """

    axes = hf.attrs["axes"]  # list of axes in the correct order

    axes_values = {
        axis: hf[f"axes/{axis}"][:] for axis in axes
    }  # dictionary of axis grid points
    # Get the properties of the grid including the dimensions etc.
    return axes, *get_grid_props_cloudy(axes, axes_values, verbose=verbose)


def check_cloudy_runs(
    grid_name,
    grid_dir,
    cloudy_dir,
    replace=False,
    files_to_check=["cont"],
):
    """
    Check that all the cloudy runs have run properly

    Arguments:
        grid_name : str
            Name of the grid
        grid_dir : str
            Parent directory for the grids.
        cloudy_dir : str
            Parent directory for the cloudy runs.
        replace : boolean
            If a run has failed simply replace the model with the previous one.
    """

    # open the new grid
    with h5py.File(f"{grid_dir}/{grid_name}.hdf5", "r") as hf:
        # Get the properties of the grid including the dimensions etc.
        (
            axes,
            n_axes,
            shape,
            n_models,
            mesh,
            model_list,
            index_list,
        ) = get_grid_properties_hf(hf)

        # list of failed models
        failed_list = []
        for i, grid_params_ in enumerate(model_list):
            infile = f"{cloudy_dir}/{grid_name}/{i+1}"
            failed = False

            # check if files exist
            for ext in files_to_check:
                if not os.path.isfile(infile + "." + ext):
                    failed = True

            # if they exist also check they have size >0
            if not failed:
                for ext in files_to_check:
                    if os.path.getsize(infile + "." + ext) < 1000:
                        failed = True

            if failed:
                print(i + 1, model_list[i])
                failed_list.append(i + 1)

                """
                If replace is specified, instead replace the grid point with
                the previous one. NOTE: this should be a last resort if the
                cloudy runs of a small number of grid points are consistently
                failing.
                """
                if replace:
                    for ext in files_to_check:
                        shutil.copyfile(
                            f"{cloudy_dir}/{grid_name}/{i}.{ext}",
                            infile + f".{ext}",
                        )

        # If the files have been replace set the failed list to empty so the
        # rest of the code can run.
        if replace:
            failed_list = []
        return failed_list


def add_spectra(grid_name, grid_dir, cloudy_dir):
    """
    Open cloudy spectra and add them to the grid.

    Arguments:
        out_grid (GridFile)
            The grid file to add the spectra to.
        cloudy_dir (str)
            Parent directory for the cloudy runs.
    """

    # The cloudy spectra to save (others can be generated later)
    spec_names = ["incident", "transmitted", "nebular", "linecont"]

    # open the new grid
    with h5py.File(f"{grid_dir}/{grid_name}.hdf5", "a") as hf:
        # Get the properties of the grid including the dimensions etc.
        (
            axes,
            n_axes,
            shape,
            n_models,
            mesh,
            model_list,
            index_list,
        ) = get_grid_properties_hf(hf)

        # Determine the cloudy version...
        cloudy_version = hf.attrs["cloudy_version"]

        # ... and use to select the correct module.
        if cloudy_version.split(".")[0] == "c23":
            cloudy = cloudy23
        elif cloudy_version.split(".")[0] == "c17":
            cloudy = cloudy17

        # Read first spectra from the first grid point to get length and
        # wavelength grid.
        lam = cloudy.read_wavelength(f"{cloudy_dir}/{grid_name}/1")

        # create a group holding the spectra in the grid file
        spectra = hf.create_group("spectra")

        # save list of spectra as attribute
        spectra.attrs["spec_names"] = spec_names

        # save the wavelength
        spectra["wavelength"] = lam

        # number of wavelength points
        nlam = len(lam)

        # make spectral grids and set them to zero
        for spec_name in spec_names:
            spectra[spec_name] = np.zeros((*shape, nlam))

        # array for holding the normalisation which is calculated below and
        # used by lines
        spectra["normalisation"] = np.ones(shape)

        for i, indices in enumerate(index_list):
            indices = tuple(indices)

            # define the infile
            infile = f"{cloudy_dir}/{grid_name}/{i+1}"

            # read the continuum file containing the spectra
            spec_dict = cloudy.read_continuum(infile, return_dict=True)

            # Calculate the specific ionising photon luminosity and use this to
            # renormalise the spectrum.
            if "log10_specific_ionising_luminosity/HI" in hf:
                # create sed object
                sed = Sed(
                    lam=lam * Angstrom,
                    lnu=spec_dict["incident"] * erg / s / Hz,
                )

                # calculate Q
                ionising_photon_production_rate = (
                    sed.calculate_ionising_photon_production_rate(
                        ionisation_energy=13.6 * eV, limit=100
                    )
                )

                # calculate normalisation
                normalisation = hf["log10_specific_ionising_luminosity/HI"][
                    indices
                ] - np.log10(ionising_photon_production_rate)

                # save normalisation for later use (rescaling lines)
                spectra["normalisation"][indices] = 10**normalisation

            # save the normalised spectrum to the correct grid point
            for spec_name in spec_names:
                spectra[spec_name][indices] = (
                    spec_dict[spec_name] * spectra["normalisation"][indices]
                )


def add_lines(
    grid_name,
    grid_dir,
    cloudy_dir,
    line_type="linelist",
    lines_to_include=False,
    include_spectra=True,
):
    """
    Open cloudy lines and add them to the HDF5 grid

    Arguments:
        grid_name : str
            Name of the grid
        grid_dir : str
            Parent directory for the grids.
        cloudy_dir : str
            Parent directory for the cloudy runs.
        line_type : str
            The type of line file to use (linelist, lines)
        dlog10_specific_ionising_lum
            The difference between the original and cloudy
            log10_specific_ionising_lum used for rescaling the cloudy spectra
    """

    # open the new grid
    with h5py.File(f"{grid_dir}/{grid_name}.hdf5", "a") as hf:
        # Get the properties of the grid including the dimensions etc.
        (
            axes,
            n_axes,
            shape,
            n_models,
            mesh,
            model_list,
            index_list,
        ) = get_grid_properties_hf(hf)

        # Determine the cloudy version...
        cloudy_version = hf.attrs["cloudy_version"]

        # ... and use to select the correct module.
        if cloudy_version.split(".")[0] == "c23":
            cloudy = cloudy23
        elif cloudy_version.split(".")[0] == "c17":
            cloudy = cloudy17

        # delete lines group if it already exists
        if "lines" in hf:
            del hf["lines"]

        # define spectra
        if include_spectra:
            spectra = hf["spectra"]
            normalisation = spectra["normalisation"][:]
            lam = spectra["wavelength"][:]

        # create group for holding lines
        lines = hf.create_group("lines")

        if line_type == "linelist":
            infile = f"{cloudy_dir}/{grid_name}/1"
            lines_to_include, _, _ = cloudy.read_linelist(
                infile, extension="emergent_elin"
            )

        # set up output arrays
        for line_id in lines_to_include:
            lines[f"{line_id}/luminosity"] = np.zeros(shape)
            lines[f"{line_id}/transmitted_continuum"] = np.zeros(shape)
            lines[f"{line_id}/nebular_continuum"] = np.zeros(shape)
            lines[f"{line_id}/continuum"] = np.zeros(shape)

        for i, indices in enumerate(index_list):
            # convert indices array to tuple
            indices = tuple(indices)

            # define the infile
            infile = f"{cloudy_dir}/{grid_name}/{i+1}"

            # get TOTAL continuum spectra
            if include_spectra:
                nebular_continuum = (
                    spectra["nebular"][indices] - spectra["linecont"][indices]
                )
                continuum = spectra["transmitted"][indices] + nebular_continuum

            # get line quantities  <--- THIS NEEDS TO CHANGE

            if line_type == "lines":
                (
                    id,
                    blend,
                    wavelength,
                    intrinsic,
                    luminosity,
                ) = cloudy.read_lines(infile)

                # identify lines we want to keep
                s = np.nonzero(np.in1d(id, np.array(lines_to_include)))[0]

                id = id[s]
                wavelength = wavelength[s]
                luminosity = luminosity[s]

            elif line_type == "linelist":
                id, wavelength, luminosity = cloudy.read_linelist(
                    infile, extension="emergent_elin"
                )

            for id_, wavelength_, luminosity_ in zip(
                id, wavelength, luminosity
            ):
                line = lines[id_]

                # save line wavelength
                line.attrs["wavelength"] = wavelength_

                # if spectra have been calculated pull the normalisation
                if include_spectra:
                    norm = normalisation[indices]
                else:
                    norm = 1.0

                # Calculate line luminosity and save it. Uses normalisation
                # from spectra.
                line["luminosity"][indices] = luminosity_ * norm  # erg s^-1

                if include_spectra:
                    # Calculate transmitted continuum at the line wavelength
                    # and save it.
                    line["transmitted_continuum"][indices] = np.interp(
                        wavelength_, lam, spectra["transmitted"][indices]
                    )  # erg s^-1 Hz^-1

                    # calculate nebular continuum at the line wavelength and
                    # save it.
                    line["nebular_continuum"][indices] = np.interp(
                        wavelength_, lam, nebular_continuum
                    )  # erg s^-1 Hz^-1

                    # calculate total continuum at the line wavelength and
                    # save it.
                    line["continuum"][indices] = np.interp(
                        wavelength_, lam, continuum
                    )  # erg s^-1 Hz^-1


if __name__ == "__main__":
    # Set up the command line arguments
    parser = Parser(
        description="Create Synthesizer HDF5 grid files from cloudy outputs."
    )

    # Path to directory where cloudy runs are
    parser.add_argument("--cloudy_dir", type=str, required=True)

    # The name of the incident grid
    parser.add_argument("--incident_grid", type=str, required=True)

    # The cloudy parameters, including any grid axes
    parser.add_argument(
        "--cloudy_params", type=str, required=False, default="c17.03-sps"
    )

    # A second cloudy parameter set which supersedes the above
    parser.add_argument(
        "--cloudy_params_addition",
        type=str,
        required=False,
    )

    # Parse arguments
    args = parser.parse_args()

    # Include spectra
    parser.add_argument(
        "--include_spectra",
        type=bool,
        default=True,
        required=False,
    )

    # Boolean flag as to whether to attempt to replace missing files
    # NOTE: this is not currently used as we should re-run cloudy or figure
    # out what went wrong when there is a failure.
    parser.add_argument("--replace", type=bool, default=False, required=False)

    # Define the line calculation method.
    parser.add_argument(
        "--line_calc_method",
        type=str,
        default="lines",
        required=False,
    )

    # Define the line calculation method.
    parser.add_argument(
        "--machine",
        type=str,
        default=None,
        required=False,
    )

    # verbosity flag
    parser.add_argument("--verbose", type=bool, required=False, default=True)

    args = parser.parse_args()

    # Unpack arguments
    grid_dir = args.grid_dir
    cloudy_dir = args.cloudy_dir
    verbose = args.verbose

    # Construct grid_name from incident grid and parameter file
    grid_name = f"{args.incident_grid}_cloudy-{args.cloudy_params}"
    incident_grid_name = args.incident_grid

    include_spectra = args.include_spectra

    # Load the cloudy parameters you are going to run
    fixed_params, grid_params = load_grid_params(args.cloudy_params)

    # If an additional parameter set is provided supersede the default
    # parameters with these and adjust the new grid name.
    if args.cloudy_params_addition:
        grid_name = (
            grid_name + f"-{args.cloudy_params_addition.split('/')[-1]}"
        )
        additional_fixed_params, additional_grid_params = load_grid_params(
            args.cloudy_params_addition
        )
        fixed_params = fixed_params | additional_fixed_params
        grid_params = grid_params | additional_grid_params

    params = fixed_params | grid_params

    # Create empty synthesizer grid
    create_empty_grid(
        grid_dir, incident_grid_name, grid_name, params, grid_params
    )

    # Check cloudy runs and potentially replace them by the nearest grid point
    # if they fail.
    failed_list = check_cloudy_runs(
        grid_name, grid_dir, cloudy_dir, replace=args.replace
    )
    print("list of failed cloudy runs:", failed_list)

    # If any runs have failed prompt to re-run.
    if len(failed_list) > 0:
        # get number of failed runs
        n_fail = len(failed_list)

        print(f"ERROR: {n_fail} cloudy runs have failed.")

        if args.machine == "apollo":
            # apollo specific command
            print("Re-run with this command:")
            print(f"qsub -t 1:{n_fail} run_grid.job")

            # replace input_names with list of failed runs
            with open(
                f"{cloudy_dir}/{grid_name}/input_names.txt", "w"
            ) as myfile:
                myfile.write("\n".join(map(str, failed_list)))

    # If no runs have failed, go ahead and add spectra and lines.
    else:
        # add spectra
        if include_spectra:
            add_spectra(grid_name, grid_dir, cloudy_dir)

        # add lines
        add_lines(
            grid_name,
            grid_dir,
            cloudy_dir,
            line_type="linelist",
            include_spectra=include_spectra,
        )
