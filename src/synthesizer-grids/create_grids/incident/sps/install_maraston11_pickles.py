"""
Download Maraston2011 and convert to HDF5 synthesizer grid.
"""
import numpy as np
import os
import argparse
from pathlib import Path
import tarfile
from unyt import erg, s, Angstrom, yr
from synthesizer.conversions import llam_to_lnu
from datetime import date
import wget
import sys

# Allow the file to use incident_utils
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from incident_utils import (
    write_data_h5py,
    write_attribute,
    add_specific_ionising_luminosity,
)  # , __tag__

# TODO: add way to automatically create /original_data/model_name and /input_data/model_name directories
# currently I'm making these manually to make the code work


def download_data(
    output_dir,
    data_url="http://www.icg.port.ac.uk/~maraston/M11/SSP_M11_Pickles.tar.gz",
):
    """
    Download Maraston+11 data
    Args:
        output_dir (string):
            directory to download and unpack data into
        data_url (string):
            URL from which to fetch the data
    Returns:
        None
    """
    filename = wget.download(
        data_url
    )  # download the original data to the working directory

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    print("filename:", filename)
    # --- untar main directory
    tar = tarfile.open(filename)
    tar.extractall(path=output_dir)
    tar.close()
    os.remove(filename)


def make_grid(model, imf, extension, output_dir):
    """Main function to convert Maraston 2011 and
    produce grids used by synthesizer
    Args:
        model (dict):
            dictionary containing model parameters
        imf (string):
            Initial mass function, can be one of Salpeter,
            Kroupa or Chabrier
        extension (string):
            String extension to use at the end of the output
            filename
        output_dir (string):
            directory where the raw Maraston+11 files are read from
    Returns:
        fname (string):
            output filename
    """

    # define output
    fname = f"{synthesizer_data_dir}/input_files/{model_name}/{model_name}{extension}_{imf}.hdf5"

    metallicities = np.array([0.02])  # array of available metallicities

    log10metallicities = np.log10(metallicities)

    metallicity_code = {0.02: "002"}  # codes for converting metallicty

    fn = f"{output_dir}/ssp_M11_Pickles{extension}.{imf_code[imf]}z{metallicity_code[metallicities[0]]}"

    ages_, _, lam_, llam_ = np.loadtxt(fn).T  # llam is in (ergs /s /AA /Msun)

    ages_Gyr = np.sort(np.array(list(set(ages_))))  # Gyr
    ages = ages_Gyr * 1e9 * yr
    log10ages = np.log10(ages)

    lam = lam_[ages_ == ages_[0]] * Angstrom

    spec = np.zeros((len(ages), len(metallicities), len(lam)))

    # at each point in spec convert the units
    for imetal, metallicity in enumerate(metallicities):
        for ia, age_Gyr in enumerate(ages_Gyr):
            print(imetal, ia, fn)
            ages_, _, lam_, llam_ = np.loadtxt(fn).T

            llam = llam_[ages_ == age_Gyr] * erg / s / Angstrom
            lnu = llam_to_lnu(lam, llam)
            spec[ia, imetal] = lnu

    # write out spectra
    write_data_h5py(fname, "spectra/wavelength", data=lam, overwrite=True)
    write_attribute(
        fname,
        "spectra/wavelength",
        "Description",
        "Wavelength of the spectra grid",
    )
    write_attribute(fname, "spectra/wavelength", "Units", "AA")

    write_data_h5py(fname, "spectra/incident", data=spec, overwrite=True)
    write_attribute(
        fname,
        "spectra/incident",
        "Description",
        "Three-dimensional spectra grid, [age, metallicity, wavelength]",
    )
    write_attribute(fname, "spectra/incident", "Units", "erg s^-1 Hz^-1")

    # write out axes
    write_attribute(fname, "/", "axes", ("log10age", "metallicity"))

    write_data_h5py(fname, "axes/log10age", data=log10ages, overwrite=True)
    write_attribute(
        fname,
        "axes/log10age",
        "Description",
        "Stellar population ages in log10 years",
    )
    write_attribute(fname, "axes/log10age", "Units", "log10(yr)")

    write_data_h5py(
        fname, "axes/metallicity", data=metallicities, overwrite=True
    )
    write_attribute(fname, "axes/metallicity", "Description", "raw abundances")
    write_attribute(
        fname, "axes/metallicity", "Units", "dimensionless [metal]"
    )

    return fname


# Lets include a way to call this script not via an entry point
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Maraston+11 download and grid creation"
    )
    parser.add_argument("-synthesizer_data_dir", type=str, required=True)
    parser.add_argument(
        "-download_data", "--download_data", type=bool, default=False
    )

    args = parser.parse_args()

    synthesizer_data_dir = args.synthesizer_data_dir

    model_name = "maraston11_pickles"

    output_dir = f"{synthesizer_data_dir}/original_data/{model_name}"  # the location to untar the original data
    imf_code = {"salpeter": "ss", "kroupa": "kr", "chabrier": "cha"}

    model = {
        "sps_name": "maraston",
        "sps_version": False,
        "alpha": False,
        "date": str(date.today()),
    }  #'synthesizer-grids_tag': __tag__,

    for extension in [
        "",
        "_nearIRextended",
        "_UVtheoretical",
        "_UVtheoretical_nearIRextended",
    ]:
        if (extension == "") or (extension == "_nearIRextended"):
            imfs = ["salpeter", "kroupa", "chabrier"]

        if (extension == "_UVtheoretical") or (
            extension == "_UVtheoretical_nearIRextended"
        ):
            imfs = ["salpeter"]

        for imf in imfs:
            print(extension)
            # if args.download_data:
            download_data(output_dir)

            fname = make_grid(
                model, imf, extension, output_dir
            )  # makes the grid and returns the name

            add_specific_ionising_luminosity(fname)
