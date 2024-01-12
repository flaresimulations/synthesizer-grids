"""
Download Maraston2011 and convert to HDF5 synthesizer grid.
"""
import numpy as np
import os
import argparse
from unyt import erg, s, Angstrom, yr
from synthesizer.conversions import llam_to_lnu
from datetime import date
import sys

# Allow the file to use incident_utils
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from incident_utils import (
    write_data_h5py,
    write_attribute,
    add_log10_specific_ionising_lum,
)  # , __tag__


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
    fname = f"{synthesizer_data_dir}/{model_name}/{model_name}_{imf}.hdf5"

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

    # at each point in spec convert the units
    for iZ, metallicity in enumerate(metallicities):
        for ia, age_Gyr in enumerate(ages_Gyr):
            print(iZ, ia, fn)
            ages_, _, lam_, llam_ = np.loadtxt(fn).T

            llam = llam_[ages_ == age_Gyr] * erg / s / Angstrom
            lnu = llam_to_lnu(lam, llam)
            spec[ia, iZ] = lnu

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
    write_attribute(fname, "axes/metallicity", "Units", "dimensionless [Z]")

    return fname


# Lets include a way to call this script not via an entry point
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Maraston+13 download and grid creation"
    )
    parser.add_argument("-synthesizer_data_dir", type=str, required=True)
    parser.add_argument(
        "-download_data", "--download_data", type=bool, default=False
    )

    args = parser.parse_args()

    synthesizer_data_dir = args.synthesizer_data_dir

    model_name = "maraston13"

    output_dir = f"{synthesizer_data_dir}/original_data/{model_name}"  # the location to untar the original data
    imfs = ["salpeter", "kroupa"]
    imf_code = {"salpeter": "ss", "kroupa": "kr"}

    model = {
        "sps_name": "maraston",
        "sps_version": False,
        "alpha": False,
        "date": str(date.today()),
    }  #'synthesizer-grids_tag': __tag__,

    for imf in imfs:
        fname = make_grid(
            model, imf, output_dir
        )  # makes the grid and returns the name

        add_log10_specific_ionising_lum(fname)
