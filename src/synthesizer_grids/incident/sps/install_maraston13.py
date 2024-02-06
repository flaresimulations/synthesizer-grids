"""
Download the Maraston2013 SPS model and convert to HDF5 synthesizer grid.
"""
import numpy as np
from unyt import erg, s, Angstrom, yr
from synthesizer.conversions import llam_to_lnu
from synthesizer_grids.parser import Parser
from synthesizer_grids.grid_io import GridFile


def make_grid(model, imf, output_dir, grid_dir):
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
        grid_dir (string):
            directory where the grids are created.
    Returns:
        fname (string):
            output filename
    """

    # define output
    out_filename = f"{grid_dir}/{model_name}_{imf}.hdf5"

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
        model=model,
        axes={"log10age": log10ages, "metallicity": metallicities},
        wavelength=lam,
        spectra={"incident": spec},  # check this unit
        alt_axes=("log10ages", "metallicities"),
    )

    # Include the specific ionising photon luminosity
    out_grid.add_specific_ionising_lum()


if __name__ == "__main__":
    # Set up the command line arguments
    parser = Parser(description="Maraston+13 download and grid creation")
    args = parser.parse_args()

    # Unpack the arguments
    args = parser.parse_args()

    # the directory to store downloaded input files
    input_dir = args.input_dir

    # the directory to store the grid
    grid_dir = args.grid_dir

    # Define the model metadata
    model_name = "maraston13"
    imfs = ["salpeter", "kroupa"]
    imf_code = {"salpeter": "ss", "kroupa": "kr"}
    model = {
        "sps_name": "maraston",
        "sps_version": False,
        "alpha": False,
    }

    # The location to untar the original data
    output_dir = f"{input_dir}/{model_name}"

    for imf in imfs:
        make_grid(model, imf, output_dir, grid_dir)
