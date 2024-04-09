"""
This reads in a cloudy grid of models

"""

import os
import shutil
from unyt import eV
import argparse
import numpy as np
import h5py

# synthesizer modules
from synthesizer.photoionisation import cloudy17, cloudy23
from synthesizer.sed import Sed

# local modules
from utils import get_grid_properties


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
    return axes, *get_grid_properties(axes, axes_values, verbose=verbose)


def check_cloudy_runs(
    grid_name,
    grid_dir,
    cloudy_dir,
    replace=False,
    files_to_check=["cont", "elin"],
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
                            infile + ".lines",
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
        grid_name (str)
            Name of the grid.
        grid_dir (str)
            Parent directory for the grids.
        cloudy_dir (str)
            Parent directory for the cloudy runs.
    """

    #  the cloudy spectra to save (others can be generated later)
    spec_names = ["incident",
                  "transmitted",
                  "nebular",
                  "linecont"]

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
        cloudy_version = hf.attrs['cloudy_version']

        # ... and use to select the correct module.
        if cloudy_version.split('.')[0] == 'c23':
            cloudy = cloudy23
        elif cloudy_version.split('.')[0] == 'c17':
            cloudy = cloudy17

        # Read first spectra from the first grid point to get length and
        # wavelength grid.
        lam = cloudy.read_wavelength(
            f"{cloudy_dir}/{grid_name}/1"
        )

        if "spectra" in hf:
            del hf["spectra"]

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
                sed = Sed(lam=lam, lnu=spec_dict["incident"])

                # calculate Q
                ionising_photon_production_rate = (
                    sed.calculate_ionising_photon_production_rate(
                        ionisation_energy=13.6*eV,
                        limit=100))

                # calculate normalisation
                normalisation = (
                    hf["log10_specific_ionising_luminosity/HI"][indices]
                    - np.log10(ionising_photon_production_rate))

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
        cloudy_version = hf.attrs['cloudy_version']

        # ... and use to select the correct module.
        if cloudy_version.split('.')[0] == 'c23':
            cloudy = cloudy23
        elif cloudy_version.split('.')[0] == 'c17':
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
            lines_to_include, _, _ = cloudy.read_linelist(infile)

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
                id, wavelength, luminosity = cloudy.read_linelist(infile)

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
    parser = argparse.ArgumentParser(
        description=("Create synthesizer HDF5 grid " "for a given grid.")
    )

    # path to grid directory (i.e. where incident and new grids are stored)
    parser.add_argument("-grid_dir", type=str, required=True)

    # path to directory where cloudy runs are
    parser.add_argument("-cloudy_dir", type=str, required=True)

    # the name of the incident grid
    parser.add_argument("-incident_grid",
                        type=str,
                        required=True)

    # the cloudy parameters, including any grid axes
    parser.add_argument("-cloudy_params",
                        type=str,
                        required=False,
                        default="c17.03-sps")

    # include spectra
    parser.add_argument(
        "-include_spectra",
        "--include_spectra",
        type=bool,
        default=True,
        required=False,
    )

    # boolean flag as to whether to attempt to replace missing files
    # NOTE: this is not currently used as we should re-run cloudy or figure
    # out what went wrong when there is a failure.
    parser.add_argument(
        "-replace", "--replace", type=bool, default=False, required=False
    )

    # Define the line calculation method.
    parser.add_argument(
        "-line_calc_method",
        "--line_calc_method",
        type=str,
        default="lines",
        required=False,
    )

    # Define the line calculation method.
    parser.add_argument(
        "-machine",
        "--machine",
        type=str,
        default=None,
        required=False,
    )

    args = parser.parse_args()

    # get arguments
    grid_dir = args.grid_dir
    cloudy_dir = args.cloudy_dir

    # construct grid_name from incident grid and parameter file
    grid_name = f"{args.incident_grid}_cloudy-{args.cloudy_params}"

    include_spectra = args.include_spectra

    # Check cloudy runs and potentially replace them by the nearest grid point
    # if they fail.
    failed_list = check_cloudy_runs(
        grid_name,
        grid_dir,
        cloudy_dir,
        replace=args.replace
    )
    print('list of failed cloudy runs:', failed_list)

    # If any runs have failed prompt to re-run.
    if len(failed_list) > 0:

        # get number of failed runs
        n_fail = len(failed_list)

        print(f"ERROR: {n_fail} cloudy runs have failed.")

        if args.machine == 'apollo':
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
            add_spectra(grid_name,
                        grid_dir,
                        cloudy_dir)

        # add lines
        add_lines(
            grid_name,
            grid_dir,
            cloudy_dir,
            line_type="linelist",
            include_spectra=include_spectra,
        )
