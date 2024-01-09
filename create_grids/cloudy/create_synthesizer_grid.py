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
from synthesizer.photoionisation import cloudy17 as cloudy
from synthesizer.sed import Sed

# local modules
from utils import get_grid_properties


def get_grid_properties_hf(hf, verbose=True):
    """
    A wrapper over get_grid_properties to get the grid properties for a HDF5 grid.
    """

    axes = hf.attrs["axes"]  # list of axes in the correct order
    axes_values = {
        axis: hf[f"axes/{axis}"][:] for axis in axes
    }  # dictionary of axis grid points

    # Get the properties of the grid including the dimensions etc.
    return axes, *get_grid_properties(axes, axes_values, verbose=verbose)


def check_cloudy_runs(
    grid_name,
    synthesizer_data_dir,
    replace=False,
    files_to_check=["cont", "elin"],
):
    """
    Check that all the cloudy runs have run properly

    Parameters
    ----------
    grid_name : str
        Name of the grid
    synthesizer_data_dir : str
        Directory where synthesizer data is kept.
    replace : boolean
        If a run has failed simply replace the model with the previous one
    """

    # open the new grid
    with h5py.File(f"{synthesizer_data_dir}/{grid_name}.hdf5", "r") as hf:
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
            infile = f"{synthesizer_data_dir}/cloudy/{grid_name}/{i+1}"
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
                If replace is specified, instead replace the grid point with the previous one.
                NOTE: this should be a last resort if the cloudy runs of a small number of grid points are consistently failing.
                """
                if replace:
                    for ext in files_to_check:
                        shutil.copyfile(
                            f"{synthesizer_data_dir}/cloudy/{grid_name}/{i}.{ext}",
                            infile + ".lines",
                        )

        # if the files have been replace set the failed list to empty so the rest of the code can run
        if replace:
            failed_list = []

        return failed_list


def add_spectra(grid_name, synthesizer_data_dir):
    """
    Open cloudy spectra and add them to the grid

    Parameters
    ----------
    grid_name : str
        Name of the grid
    synthesizer_data_dir : str
        Directory where synthesizer data is kept.
    """

    # spec_names = ['incident','transmitted','nebular','nebular_continuum','total','linecont']
    #  the cloudy spectra to save (others can be generated later)
    spec_names = ["incident", "transmitted", "nebular", "linecont"]

    # open the new grid
    with h5py.File(f"{synthesizer_data_dir}/{grid_name}.hdf5", "a") as hf:
        # Get the properties of the grid including the dimensions etc.
        (
            axes,
            n_axes,
            shape,
            n_models,
            mesh,
            model_list,
            index_list,
        ) = get_grid_properties_hf(hf, verbose=False)

        # read first spectra from the first grid point to get length and wavelength grid
        lam = cloudy.read_wavelength(
            f"{synthesizer_data_dir}/cloudy/{grid_name}/1"
        )

        if "spectra" in hf:
            del hf["spectra"]

        spectra = hf.create_group(
            "spectra"
        )  # create a group holding the spectra in the grid file
        spectra.attrs[
            "spec_names"
        ] = spec_names  # save list of spectra as attribute

        spectra["wavelength"] = lam  # save the wavelength
        nu = 3e8 / (lam * 1e-10)

        nlam = len(lam)  # number of wavelength points

        # make spectral grids and set them to zero
        for spec_name in spec_names:
            spectra[spec_name] = np.zeros((*shape, nlam))

        # array for holding the normalisation which is calculated below and used by lines
        spectra["normalisation"] = np.ones(shape)

        for i, indices in enumerate(index_list):
            indices = tuple(indices)

            # define the infile
            infile = f"{synthesizer_data_dir}/cloudy/{grid_name}/{i+1}"

            # read the continuum file containing the spectra
            spec_dict = cloudy.read_continuum(infile, return_dict=True)

            # if hf['log10Q/HI'] already exists (only true for non-cloudy
            # models) renormalise the spectrum.

            # need to add sed object and update calculate_Q

            if "log10Q/HI" in hf:
                # create sed object
                sed = Sed(lam=lam, lnu=spec_dict["incident"])

                # calculate Q
                Q = sed.calculate_ionising_photon_production_rate(
                    ionisation_energy=13.6 * eV, limit=100
                )

                # Q = calculate_Q(lam,
                #                spec_dict['incident'],
                #                ionisation_energy=13.6 * eV)

                # calculate normalisation
                normalisation = hf["log10Q/HI"][indices] - np.log10(Q)

                # save normalisation for later use (rescaling lines)
                spectra["normalisation"][indices] = 10**normalisation

            # save the normalised spectrum to the correct grid point
            for spec_name in spec_names:
                spectra[spec_name][indices] = (
                    spec_dict[spec_name] * spectra["normalisation"][indices]
                )


def add_lines(
    grid_name,
    synthesizer_data_dir,
    line_type="linelist",
    lines_to_include=False,
    include_spectra=True,
):
    """
    Open cloudy lines and add them to the HDF5 grid

    Parameters
    ----------
    grid_name: str
        Name of the grid.
    synthesizer_data_dir: str
        Directory where synthesizer data is kept.
    line_type : str
        The type of line file to use (linelist, lines)
    dlog10Q
        The difference between the original and cloudy log10Q used for rescaling the cloudy spectra
    """

    # open the new grid
    with h5py.File(f"{synthesizer_data_dir}/{grid_name}.hdf5", "a") as hf:
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
        # lines.attrs['lines'] = list(lines_to_include)  # save list of spectra as attribute

        if line_type == "linelist":
            infile = f"{synthesizer_data_dir}/cloudy/{grid_name}/1"
            lines_to_include, _, _ = cloudy.read_linelist(infile)
            # print(lines_to_include)

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
            infile = f"{synthesizer_data_dir}/cloudy/{grid_name}/{i+1}"

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

                # calculate line luminosity and save it. Uses normalisation from spectra.
                line["luminosity"][indices] = luminosity_ * norm  # erg s^-1

                if include_spectra:
                    # calculate transmitted continuum at the line wavelength and save it.
                    line["transmitted_continuum"][indices] = np.interp(
                        wavelength_, lam, spectra["transmitted"][indices]
                    )  # erg s^-1 Hz^-1

                    # calculate nebular continuum at the line wavelength and save it.
                    line["nebular_continuum"][indices] = np.interp(
                        wavelength_, lam, nebular_continuum
                    )  # erg s^-1 Hz^-1

                    # calculate total continuum at the line wavelength and save it.
                    line["continuum"][indices] = np.interp(
                        wavelength_, lam, continuum
                    )  # erg s^-1 Hz^-1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=("Create synthesizer HDF5 grid " "for a given grid.")
    )

    parser.add_argument(
        "-synthesizer_data_dir", type=str, required=True
    )  # path to synthesizer_data_dir
    parser.add_argument("-grid_name", "--grid_name", type=str, required=True)
    parser.add_argument(
        "-include_spectra",
        "--include_spectra",
        type=bool,
        default=True,
        required=False,
    )
    parser.add_argument(
        "-replace", "--replace", type=bool, default=False, required=False
    )
    parser.add_argument(
        "-line_calc_method",
        "--line_calc_method",
        type=str,
        default="lines",
        required=False,
    )

    args = parser.parse_args()

    # get arguments
    synthesizer_data_dir = args.synthesizer_data_dir
    grid_name = args.grid_name
    include_spectra = args.include_spectra

    # check cloudy runs and potentially replace them by the nearest grid point if they fail.
    failed_list = check_cloudy_runs(
        grid_name, synthesizer_data_dir, replace=args.replace
    )

    print(failed_list)

    # if failed prompt to re-run
    if len(failed_list) > 0:
        print(
            f"ERROR: {len(failed_list)} cloudy runs have failed. You should re-run these with command:"
        )
        print(f"  qsub -t 1:{len(failed_list)}  run_grid.job")

        # replace input_names with list of failed runs
        with open(
            f"{synthesizer_data_dir}/cloudy/{grid_name}/input_names.txt", "w"
        ) as myfile:
            myfile.write("\n".join(map(str, failed_list)))

    # if not failed, go ahead and add spectra and lines
    else:
        print("- passed checks")

        # add spectra
        if include_spectra:
            add_spectra(grid_name, synthesizer_data_dir)
            print("- spectra added")

        # add lines
        add_lines(
            grid_name,
            synthesizer_data_dir,
            line_type="linelist",
            include_spectra=include_spectra,
        )
