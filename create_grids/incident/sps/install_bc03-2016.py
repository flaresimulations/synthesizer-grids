"""
Download BC03 and convert to HDF5 synthesizer grid.

TODO: data download does not currently work, directories need to be updated
"""
import os
import sys
import argparse
import numpy as np
import re
import wget
import tarfile
import glob
import gzip
import shutil
from datetime import date
from unyt import angstrom, erg, s, Hz

sys.path.insert(1, os.path.dirname(os.path.abspath(sys.argv[0])) + "/../../")
from grid_io import GridFile
from utils import (
    get_model_filename,
)
from synthesizer._version import __version__


def download_data(variant):
    url = (
        "http://www.bruzual.org/bc03/Updated_version_2016/"
        f"BC03_{variant.lower()}_chabrier.tgz"
    )

    filename = wget.download(url)
    return filename


def untar_data(synthesizer_data_dir):
    input_dir = f"{synthesizer_data_dir}/input_files/"
    fn = "bc03.models.padova_2000_chabrier_imf.tar.gz"

    # --- untar main directory
    tar = tarfile.open(fn)
    tar.extractall(path=input_dir)
    tar.close()
    os.remove(fn)

    # --- unzip the individual files that need reading
    model_dir = (
        f"{synthesizer_data_dir}/input_files/"
        "bc03/models/Padova2000/chabrier"
    )
    files = glob.glob(f"{model_dir}/bc2003_hr_m*_chab_ssp.ised_ASCII.gz")

    for file in files:
        with gzip.open(file, "rb") as f_in:
            with open(".".join(file.split(".")[:-1]), "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)


def readBC03Array(file, lastLineFloat=None):
    """Read a record from bc03 ascii file. The record starts with the
       number of elements N and is followed by N numbers. The record may
       or may not start within a line, i.e. a line need not necessarily
       start with a record.
    Parameters:
    ----------------------------------------------------------------------
    file: handle on open bc03 ascii file
    lastLineFloat: still open line from last line read, in case of a
                   record starting mid-line.
    Returns array, lastLine, where:
    ----------------------------------------------------------------------
    array = The array values read from the file
    lastLine = The remainder of the last line read (in floating format),
               for continued reading of the file
    """

    if lastLineFloat is None or len(lastLineFloat) == 0:
        # Nothing in last line, so read next line
        line = file.readline()
        lineStr = line.split()
        lastLineFloat = [float(x) for x in lineStr]
    # Read array 'header' (i.e. number of elements)
    arrayCount = int(lastLineFloat[0])  # Length of returned array
    array = np.empty(arrayCount)  # Initialise the array
    lastLineFloat = lastLineFloat[1 : len(lastLineFloat)]
    iA = 0  # Running array index
    while True:  # Read numbers until array is full
        for iL in range(0, len(lastLineFloat)):  # Loop numbers in line
            array[iA] = lastLineFloat[iL]
            iA = iA + 1
            if iA >= arrayCount:  # Array is full so return
                return array, lastLineFloat[iL + 1 :]
        line = file.readline()  # Went through the line so get the next one
        lineStr = line.split()
        lastLineFloat = [float(x) for x in lineStr]


def convertBC03(files=None):
    """Convert BC03 outputs
    Parameters (user will be prompted for those if not present):
    ----------------------------------------------------------------------
    files: list of each BC03 SED ascii file, typically named
           bc2003_xr_mxx_xxxx_ssp.ised_ASCII
    """

    # Prompt user for files if not provided--------------------
    if files is None:
        print(
            "Please write the model to read",
        )
        files = []
        while True:
            filename = input("filename >")
            if filename == "":
                break
            files.append(filename)
        print("ok checking now")
        if not len(files):
            print("No BC03 files given, nothing do to")
            return

    # Initialise ---------------------------------------------------------
    ageBins = None
    lambdaBins = None
    metalBins = [None] * len(files)
    seds = np.array([[[None]]])

    print("Reading BC03 files and converting...")
    # Loop SED tables for different metallicities
    for iFile, fileName in enumerate(files):
        print("Converting file ", fileName)
        file = open(fileName, "r")
        # file = gzip.open(f'{fileName}.gz', 'rb')

        ages, lastLine = readBC03Array(file)  # Read age bins
        nAge = len(ages)
        print("Number of ages: %s" % nAge)
        if ageBins is None:
            ageBins = ages
            seds.resize(
                (seds.shape[0], len(ageBins), seds.shape[2]), refcheck=False
            )
        if not np.array_equal(ages, ageBins):  # check for consistency
            print("Age bins are not identical everywhere!!!")
            print("CANCELLING CONVERSION!!!")
            return
        # Read four (five ?) useless lines
        line = file.readline()
        line = file.readline()
        line = file.readline()
        line = file.readline()
        line = file.readline()
        # These last three lines are identical and contain the metallicity
        (jmetal,) = re.search("metal=([0-9]+\.?[0-9]*)", line).groups()
        metalBins[iFile] = eval(jmetal)
        seds.resize(
            (len(metalBins), seds.shape[1], seds.shape[2]), refcheck=False
        )
        # Read wavelength bins
        lambdas, lastLine = readBC03Array(file, lastLineFloat=lastLine)
        if lambdaBins is None:  # Write wavelengths to sed file
            lambdaBins = lambdas
            seds.resize(
                (seds.shape[0], seds.shape[1], len(lambdaBins)), refcheck=False
            )
        if not np.array_equal(lambdas, lambdaBins):  # check for consistency
            print("Wavelength bins are not identical everywhere!!!")
            print("CANCELLING CONVERSION!!!")
            return
        # Read luminosities
        for iAge in range(0, nAge):
            lums, lastLine = readBC03Array(file, lastLineFloat=lastLine)
            if len(lums) != len(lambdaBins):
                print("Inconsistent number of wavelength bins in BC03")
                print("STOPPING!!")
                return
            # Read useless array
            tmp, lastLine = readBC03Array(file, lastLineFloat=lastLine)
            seds[iFile, iAge] = lums
            progress = (iAge + 1) / nAge
            sys.stdout.write(
                "\rProgress: [{0:50s}] {1:.1f}%".format(
                    "#" * int(progress * 50), progress * 100
                )
            )
        print(" ")
        lastLine = None

    return (
        np.array(seds, dtype=np.float64),
        np.array(metalBins, dtype=np.float64),
        np.array(ageBins, dtype=np.float64),
        np.array(lambdaBins, dtype=np.float64),
    )


def make_grid(variant, synthesizer_data_dir, out_filename):
    """Main function to convert BC03 grids and
    produce grids used by synthesizer"""

    # Define base path
    if variant == "BaSeL":
        variant_dir = "BaSeL3.1_Atlas"
        variant_code = "lr_BaSeL"
    if variant == "Miles":
        variant_dir = "Miles_Atlas"
        variant_code = "hr_xmiless"
    if variant == "Stelib":
        variant_dir = "Stelib_Atlas"
        variant_code = "hr_stelib"

    basepath = (
        f"{synthesizer_data_dir}/input_files/"
        f"bc03-2016/{variant_dir}/Chabrier_IMF/"
    )

    # Define files
    files = [
        f"bc2003_{variant_code}_m22_chab_ssp.ised_ASCII",
        f"bc2003_{variant_code}_m32_chab_ssp.ised_ASCII",
        f"bc2003_{variant_code}_m42_chab_ssp.ised_ASCII",
        f"bc2003_{variant_code}_m52_chab_ssp.ised_ASCII",
        f"bc2003_{variant_code}_m62_chab_ssp.ised_ASCII",
        f"bc2003_{variant_code}_m72_chab_ssp.ised_ASCII",
        f"bc2003_{variant_code}_m82_chab_ssp.ised_ASCII",
    ]

    out = convertBC03([basepath + s for s in files])

    metallicities = out[1]
    log10metallicities = np.log10(metallicities)

    ages = out[2]
    ages[0] = 1e5
    log10ages = np.log10(ages)

    lam = out[3]
    nu = 3e8 / (lam * 1e-10)

    spec = out[0]

    spec = np.swapaxes(spec, 0, 1)  # make (age, metallicity, wavelength)

    spec *= 3.826e33  # erg s^-1 AA^-1 Msol^-1
    spec *= lam / nu  # erg s^-1 Hz^-1 Msol^-1

    na = len(ages)
    nmetal = len(metallicities)

    # Create the GridFile ready to take outputs
    out_grid = GridFile(out_filename, mode="a", overwrite=True)

    # Write everything out thats common to all models
    out_grid.write_grid_common(
        model,
        axes={"log10age": log10ages, "metallicity": metallicities},
        wavelength=lam * angstrom,
        spectra={"incident": spec * erg / s / Hz},
        alt_axes=("log10ages", "metallicities"),
    )

    # Include the specific ionising photon luminosity
    out_grid.add_specific_ionising_lum()

    out_grid.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="BC03-2016 download and grid creation"
    )
    parser.add_argument(
        "-synthesizer_data_dir",
        "--synthesizer_data_dir",
        default="/Users/sw376/Dropbox/Research/data/synthesizer",
    )
    parser.add_argument(
        "-download_data", "--download_data", type=bool, default=False
    )

    args = parser.parse_args()

    synthesizer_data_dir = args.synthesizer_data_dir

    grid_dir = f"{synthesizer_data_dir}/grids"

    # download data
    if args.download_data:
        download_data()
        untar_data()

    default_model = {
        "sps_name": "bc03-2016",
        "sps_version": False,
        "sps_variant": False,
        "imf_type": "chabrier03",  # named IMF or bpl (broken power law)
        "imf_masses": [0.1, 100],
        "imf_slopes": False,
        "alpha": False,
        "synthesizer-grids_tag": __version__,
        "date": str(date.today()),
    }

    for variant in ["BaSeL", "Miles", "Stelib"]:  # 'BaSeL',
        model = default_model | {"sps_variant": variant}

        synthesizer_model_name = get_model_filename(model)
        print(synthesizer_model_name)

        # this is the full path to the ultimate HDF5 grid file
        out_filename = (
            f"{synthesizer_data_dir}/grids/dev/{synthesizer_model_name}.hdf5"
        )

        make_grid(variant, synthesizer_data_dir, out_filename)
