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

from synthesizer_grids.utilities import GridFile
from sythesizer_grids.utilities import Parser
from incident_utils import (
    write_data_h5py,
    write_attribute,
    add_log10_specific_ionising_lum,
)  # , __tag__


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
    out_filename = (
        f"{synthesizer_data_dir}/grids/{model_name}{extension}_{imf}.hdf5"
    )

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

    # Create the GridFile ready to take outputs
    out_grid = GridFile(out_filename, mode="w", overwrite=True)

    # at each point in spec convert the units
    for imetal, metallicity in enumerate(metallicities):
        for ia, age_Gyr in enumerate(ages_Gyr):
            print(imetal, ia, fn)
            ages_, _, lam_, llam_ = np.loadtxt(fn).T

            llam = llam_[ages_ == age_Gyr] * erg / s / Angstrom
            lnu = llam_to_lnu(lam, llam)
            spec[ia, imetal] = lnu

        # Write everything out thats common to all models
    out_grid.write_grid_common(
        model,
        axes={"log10age": log10ages, "metallicity": metallicities},
        wavelength=lam,
        spectra={"incident": spec},
        alt_axes=("log10ages", "metallicities"),
    )

    # Include the specific ionising photon luminosity
    out_grid.add_specific_ionising_lum()


# Lets include a way to call this script not via an entry point
if __name__ == "__main__":
    # Set up the command line arguments
    parser = Parser(description="Maraston+11 download and grid creation")
    args = parser.parse_args()

    # Unpack the arguments
    synthesizer_data_dir = args.synthesizer_data_dir
    grid_dir = f"{synthesizer_data_dir}/grids"

    # Define the model metadata
    model_name = "maraston11_pickles"
    imf_code = {"salpeter": "ss", "kroupa": "kr", "chabrier": "cha"}
    model = {
        "sps_name": "maraston",
        "sps_version": False,
        "alpha": False,
        "date": str(date.today()),
    }

    # The location to untar the original data
    output_dir = f"{synthesizer_data_dir}/original_data/{model_name}"

    # Download the data if necessary
    if args.download:
        download_data(output_dir)

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

            make_grid(
                model, imf, extension, output_dir
            )  # makes the grid and returns the name
