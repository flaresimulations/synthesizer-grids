"""
Download BC03 and convert to HDF5 synthesizer grid.
"""

import numpy as np
import os
import wget
import argparse
from pathlib import Path
import tarfile
from synthesizer.utils import flam_to_fnu
from synthesizer.sed import calc_specific_ionising_luminosity
from datetime import date
from utils import (
    __tag__,
    write_data_h5py,
    write_attribute,
    add_specific_ionising_luminosity,
    get_model_filename,
)


def download_data(input_dir):
    """
    TODO: These could be replaced by our own mirror
    """

    filename = wget.download(
        original_data_url[imf]
    )  # download the original data to the working directory

    Path(input_dir).mkdir(parents=True, exist_ok=True)

    # --- untar main directory
    tar = tarfile.open(filename)
    tar.extractall(path=input_dir)
    tar.close()
    os.remove(filename)

    # # --- unzip the individual files that need reading
    # model_dir = f'{sythesizer_data_dir}/input_files/bc03/models/Padova2000/chabrier'
    # files = glob.glob(f'{model_dir}/bc2003_hr_m*_chab_ssp.ised_ASCII.gz')
    #
    # for file in files:
    #     with gzip.open(file, 'rb') as f_in:
    #         with open('.'.join(file.split('.')[:-1]), 'wb') as f_out:
    #             shutil.copyfileobj(f_in, f_out)
    #


def make_grid(model, imf, hr_morphology):
    """Main function to convert BC03 grids and
    produce grids used by synthesizer"""

    # get synthesizer model name
    synthesizer_model_name = get_model_filename(model)

    print(synthesizer_model_name)

    # Define output
    out_filename = (
        f"{synthesizer_data_dir}/grids/dev/{synthesizer_model_name}.hdf5"
    )

    # NOTE THE LOWEST METALLICITY MODEL DOES NOT HAVE YOUNG AGES so don't use
    metallicities = np.array(
        [0.001, 0.01, 0.02, 0.04]
    )  # array of avialable metallicities

    # codes for converting metallicty
    metallicity_code = {
        0.0001: "10m4",
        0.001: "0001",
        0.01: "001",
        0.02: "002",
        0.04: "004",
        0.07: "007",
    }

    # --- open first file to get age
    fn = f"{input_dir}/sed.{imf}z{metallicity_code[metallicities[0]]}.{hr_morphology}"
    ages_, _, lam_, flam_ = np.loadtxt(fn).T

    ages_Gyr = np.sort(np.array(list(set(ages_))))  # Gyr
    ages = ages_Gyr * 1e9  # yr
    log10ages = np.log10(ages)

    lam = lam_[ages_ == ages_[0]]

    spec = np.zeros((len(ages), len(metallicities), len(lam)))

    for imetal, metallicity in enumerate(metallicities):
        for ia, age_Gyr in enumerate(ages_Gyr):
            fn = f"{input_dir}/sed.{imf}z{metallicity_code[metallicity]}.{hr_morphology}"
            print(imetal, ia, fn)
            ages_, _, lam_, flam_ = np.loadtxt(fn).T

            flam = flam_[ages_ == age_Gyr]
            fnu = flam_to_fnu(lam, flam)
            spec[ia, imetal] = fnu

    write_data_h5py(
        out_filename, "spectra/wavelength", data=lam, overwrite=True
    )
    write_attribute(
        out_filename,
        "spectra/wavelength",
        "Description",
        "Wavelength of the spectra grid",
    )
    write_attribute(out_filename, "spectra/wavelength", "Units", "AA")

    write_data_h5py(
        out_filename, "spectra/incident", data=spec, overwrite=True
    )
    write_attribute(
        out_filename,
        "spectra/incident",
        "Description",
        "Three-dimensional spectra grid, [age, metallicity, wavelength]",
    )
    write_attribute(out_filename, "spectra/incident", "Units", "erg/s/Hz")

    # write out axes
    write_attribute(out_filename, "/", "axes", ("log10age", "metallicity"))

    write_data_h5py(
        out_filename, "axes/log10age", data=log10ages, overwrite=True
    )
    write_attribute(
        out_filename,
        "axes/log10age",
        "Description",
        "Stellar population ages in log10 years",
    )
    write_attribute(out_filename, "axes/log10age", "Units", "dex(yr)")

    write_data_h5py(
        out_filename, "axes/metallicity", data=metallicities, overwrite=True
    )
    write_attribute(
        out_filename, "axes/metallicity", "Description", "raw abundances"
    )
    write_attribute(out_filename, "axes/metallicity", "Units", "dimensionless")

    return out_filename


# Lets include a way to call this script not via an entry point
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Maraston download and grid creation"
    )
    parser.add_argument(
        "--download",
        default=False,
        action="store_true",
        help=(
            "download bpass data directly in current directory "
            "and untar in sunthesizer data directory"
        ),
    )

    parser.add_argument(
        "-synthesizer_data_dir",
        "--synthesizer_data_dir",
        default="/Users/sw376/Dropbox/Research/data/synthesizer",
    )

    args = parser.parse_args()

    synthesizer_data_dir = args.synthesizer_data_dir

    model_name = "maraston"

    input_dir = f"{synthesizer_data_dir}/input_files/{model_name}"  # the location to untar the original data

    imfs = ["ss"]  # , 'kr'

    original_data_url = {}
    original_data_url[
        "ss"
    ] = "http://www.icg.port.ac.uk/~maraston/SSPn/SED/Sed_Mar05_SSP_Salpeter.tar.gz"

    model = {
        "sps_name": "maraston",
        "sps_version": False,
        "alpha": False,
        "synthesizer-grids_tag": __tag__,
        "date": str(date.today()),
    }

    for imf in imfs:
        if imf == "ss":
            model["imf_type"] = "bpl"
            model["imf_masses"] = [0.1, 100]
            model["imf_slopes"] = [2.35]

        for hr_morphology in ["rhb"]:
            model["sps_variant"] = hr_morphology
            out_filename = make_grid(model, imf, hr_morphology)

            add_specific_ionising_luminosity(out_filename)
