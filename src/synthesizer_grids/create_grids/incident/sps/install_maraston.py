"""
Download BC03 and convert to HDF5 synthesizer grid.
"""
import numpy as np
import os
import wget
from pathlib import Path
import tarfile
from synthesizer.utils import flam_to_fnu
from datetime import date

from synthesizer_grids.utilities.parser import Parser
from synthesizer_grids.utilities import GridFile
from utils import get_model_filename


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

    # Create the GridFile ready to take outputs
    out_grid = GridFile(out_filename, mode="w", overwrite=True)

    # Write everything out thats common to all models
    out_grid.write_grid_common(
        model=model,
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
    parser = Parser(description="Maraston download and grid creation")
    args = parser.parse_args()

    # Unpack the arguments
    synthesizer_data_dir = args.synthesizer_data_dir
    grid_dir = f"{synthesizer_data_dir}/grids"

    # Define the model metadata
    model_name = "maraston"
    imfs = ["ss"]  # , 'kr'
    model = {
        "sps_name": "maraston",
        "sps_version": False,
        "alpha": False,
    }

    # Define the download URL
    original_data_url = {}
    original_data_url[
        "ss"
    ] = "http://www.icg.port.ac.uk/~maraston/SSPn/SED/Sed_Mar05_SSP_Salpeter.tar.gz"

    # The location to untar the original data
    input_dir = f"{synthesizer_data_dir}/input_files/{model_name}"

    for imf in imfs:
        if imf == "ss":
            model["imf_type"] = "bpl"
            model["imf_masses"] = [0.1, 100]
            model["imf_slopes"] = [2.35]

        for hr_morphology in ["rhb"]:
            model["sps_variant"] = hr_morphology

            # Get and write the grid
            make_grid(model, imf, hr_morphology)
