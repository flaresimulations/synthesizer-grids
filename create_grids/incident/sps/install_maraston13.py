"""
Download the Maraston2013 SPS model and convert to HDF5 synthesizer grid.
"""
import numpy as np
import os
import argparse
from unyt import erg, s, Angstrom, yr
from synthesizer.conversions import llam_to_lnu
from datetime import date
import sys

# to allow access to the grid_io module:
sys.path.insert(1, os.path.dirname(os.path.abspath(sys.argv[0])) + "/../../")

from grid_io import GridFile


def make_grid(model, imf, output_dir):
    """Main function to convert Maraston 2013 and
    produce grids used by synthesizer
    Args:
        model (dict):
            dictionary containing model parameters
        imf (string):
            Initial mass function, Salpeter or
            Kroupa
        output_dir (string):
            directory where the raw Maraston+13 files are read from
    Returns:
        fname (string):
            output filename
    """

    # define output
    out_filename = f"{synthesizer_data_dir}/grids/{model_name}_{imf}.hdf5"

    metallicities = np.array(
        [0.01, 0.001, 0.02, 0.04]
    )  # array of available metallicities

    metallicity_code = {
        0.01: "001",
        0.001: "0001",
        0.02: "002",
        0.04: "004",
    }  # codes for converting metallicty

    # open first raw data file to get age
    fn = f"{output_dir}/sed_M13.{imf_code[imf]}z{metallicity_code[metallicities[0]]}"

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
        spectra={"incident": spec},  # check this unit
        alt_axes=("log10ages", "metallicities"),
    )

    # Include the specific ionising photon luminosity
    out_grid.add_specific_ionising_lum()


# Lets include a way to call this script not via an entry point
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Maraston+13 download and grid creation"
    )
    parser.add_argument("-synthesizer_data_dir", type=str, required=True)
    parser.add_argument("-download_data", "--download_data", type=bool, default=False)

    args = parser.parse_args()

    synthesizer_data_dir = args.synthesizer_data_dir

    model_name = "maraston13"

    output_dir = f"{synthesizer_data_dir}/input_files/{model_name}"  # the location to untar the original data
    imfs = ["salpeter", "kroupa"]
    imf_code = {"salpeter": "ss", "kroupa": "kr"}

    model = {
        "sps_name": "maraston",
        "sps_version": False,
        "alpha": False,
        "date": str(date.today()),
    }  #'synthesizer-grids_tag': __tag__,

    for imf in imfs:
        make_grid(model, imf, output_dir)
