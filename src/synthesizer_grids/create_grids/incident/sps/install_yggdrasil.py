"""
Download Populations 3 star spectra from Yggdrasil models and convert to HDF5
synthesizer grid.

There are 3 versions implemented:
    1. Pop III.1 - A zero-metallicity population with an extremely top-heavy
                   IMF (50-500 Msolar, Salpeter slope), using a SSP from
                   Schaerer et al. (2002, A&A, 382, 28)
    2. Pop III.2 - A zero-metallicity population with a moderately top-heavy
                   IMF (log-normal with characteristic mass M_c=10 Msolar,
                   dispersion sigma=1 Msolar and wings extending from 1-500
                   Msolar) from Raiter et al. (2010, A&A 523, 64)
    3. Pop III, Kroupa IMF - A zero-metallicity population with a normal
                             IMF (universal Kroupa 2001 IMF in the interval
                             0.1-100 Msolar), based on a rescaled SSP from
                             Schaerer et al. (2002, A&A, 382, 28)

We also just pick the instantaneous burst model with the 3 different gas
covering factor of 0 (no nebular contribution), 0.5, 1 (maximal nebular
contribution)

Warning: the nebular procesed grids here differ from the rest of the
nebular processing implementation in synthesizer, where we self consistently
run pure stellar spectra through CLOUDY. For full self consistency the
nebular grids here should not be used, but we provide anyway for reference.
"""
import numpy as np
import re
import requests
from tqdm import tqdm
from unyt import erg, s, angstrom, c, Hz

from synthesizer_grids.utilities.parser import Parser
from synthesizer_grids.utilities.grid_io import GridFile


def download_data(synthesizer_data_dir, ver, fcov):
    """
    Function access Yggdrasil spectra from website
    """
    # Define base path
    input_dir = f"{synthesizer_data_dir}/input_files/"
    filename = f"PopIII{ver}_fcov_{fcov}_SFR_inst_Spectra"
    save_path = f"{input_dir}/popIII/Yggdrasil/{filename}"

    # Define the url
    url = f"https://www.astro.uu.se/~ez/yggdrasil/YggdrasilSpectra/{filename}"

    # Call the server and get the response
    response = requests.get(url, stream=True, verify=False)

    # Sizes in bytes.
    total_size = int(response.headers.get("content-length", 0))
    block_size = 1024

    # Stream the file to disk with a nice progress bar.
    with tqdm(total=total_size, unit="B", unit_scale=True) as progress_bar:
        with open(save_path, "wb") as f:
            for chunk in response.iter_content(block_size):
                progress_bar.update(len(chunk))
                f.write(chunk)

    return save_path


def convertPOPIII(synthesizer_data_dir, ver, fcov, fileloc):
    """
    Convert POPIII outputs for Yggdrasil
    Wavelength in Angstrom
    Flux is in erg/s/AA
    """

    # Initialise ---------------------------------------------------------
    ageBins = None
    lambdaBins = None
    metalBins = np.array([0])
    seds = np.array([[[None]]])

    print("Reading POPIII files and converting...")
    # Open SED table containing different ages
    print("Converting file ", fileloc)
    data = open(fileloc, "r")
    text = data.read()

    # Get age values
    ages = re.findall(r"Age\s(.*?)\n", text)
    ageBins = np.array(
        [
            float(
                re.findall(
                    r" [+\-]?[^\w]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)",
                    ages[ii],
                )[0]
            )
            for ii in range(len(ages))
        ]
    )

    # Get the number of wavelength points: lam_num
    lam_num = np.array(re.findall(r"points:(.*?)\n", text), dtype=int)
    diff = np.diff(lam_num)
    if np.sum(diff) != 0:
        print("Age bins are not identical everywhere!!!")
        print("CANCELLING CONVERSION!!!")
        return

    seds = np.zeros((len(ageBins), len(metalBins), lam_num[0]))

    """ 
        Format of the file is 10 header lines at begining followed by
        lam_num lines of wavelength and flux, then one empty line and
        7 string lines giving the ages 
    """
    data = open(fileloc, "r")
    tmp = data.readlines()
    mass = float(re.findall(r"\d+\.\d+", tmp[0])[0])
    begin = 9
    end = begin + lam_num[0]
    for ii in range(len(ageBins)):
        this_data = tmp[begin:end]
        if ii == 0:
            lambdaBins = np.array(
                [
                    float(re.findall(r"[-+]?([0-9]*\.[0-9]+|[0-9]+)", jj)[0])
                    for jj in this_data
                ]
            )

        seds[ii, 0] = np.array(
            [
                float(re.findall(r"[-+]?([0-9]*\.[0-9]+|[0-9]+)", jj)[1])
                * (
                    10
                    ** float(
                        re.findall(r"[-+]?([0-9]*\.[0-9]+|[0-9]+)", jj)[2]
                    )
                )
                for jj in this_data
            ]
        )

        begin = end + 8
        end = begin + lam_num[0]

    return (
        np.array(seds / mass, dtype=np.float64),
        np.array(metalBins, dtype=np.float64),
        np.array(ageBins, dtype=np.float64),
        np.array(lambdaBins, dtype=np.float64),
    )


def make_grid(synthesizer_data_dir, ver, fcov):
    """Main function to convert POPIII grids and
    produce grids used by synthesizer"""

    model_name = f"yggdrasil_POPIII{ver}"
    out_filename = f"{synthesizer_data_dir}/grids/{model_name}.hdf5"

    # Define input path
    input_dir = f"{synthesizer_data_dir}/input_files/"
    filename = f"PopIII{ver}_fcov_{fcov}_SFR_inst_Spectra"
    input_path = f"{input_dir}/popIII/Yggdrasil/{filename}"

    # Get spectra and attributes
    out = convertPOPIII(synthesizer_data_dir, ver, fcov, input_path)

    metallicities = out[1]

    ages = out[2] * 1e6  # since ages are quoted in Myr
    log10ages = np.log10(ages)

    lam = out[3]

    # Converting L_lam to L_nu using
    # L_lam dlam = L_nu dnu
    # L_nu = L_lam (lam)^2 / c
    # c in units of AA/s for conversion

    light_speed = c.to(angstrom / s).value  # in AA/s
    spec = out[0]

    spec *= (lam**2) / light_speed  # now in erg s^-1 Hz^-1 Msol^-1

    # Create the grid file
    out_grid = GridFile(out_filename, mode="a", overwrite=True)

    if fcov == "0":
        # Write everything out thats common to all models
        out_grid.write_grid_common(
            axes={"log10age": log10ages, "metallicity": metallicities},
            wavelength=lam * angstrom,
            spectra={"incident": spec * erg / s / Hz},
            alt_axes=("log10ages", "metallicities"),
        )

    else:
        # Do we needa  suffix?
        if fcov == "1":
            add = ""
        else:
            add = f"_fcov_{fcov}"

        # Write everything out thats common to all models
        out_grid.write_grid_common(
            axes={"log10age": log10ages, "metallicity": metallicities},
            wavelength=lam * angstrom,
            spectra={f"nebular{add}": spec * erg / s / Hz},
            alt_axes=("log10ages", "metallicities"),
        )

    # Include the specific ionising photon luminosity
    out_grid.add_specific_ionising_lum()


if __name__ == "__main__":
    # Set up the command line arguments
    parser = Parser(description="Yggdrasil download and grid creation")
    args = parser.parse_args()

    # Unpack the arguments
    synthesizer_data_dir = args.synthesizer_data_dir
    grid_dir = f"{synthesizer_data_dir}/grids"

    # Different forms of the IMFs
    vers = np.array([".1", ".2", "_kroupa_IMF"])

    # Different gas covering fractions for nebular emission model
    fcovs = np.array(["0", "0.5", "1"])

    for ver in vers:
        for fcov in fcovs:
            # Download the data if necessary
            if args.download:
                download_data(synthesizer_data_dir, ver, fcov)

            make_grid(synthesizer_data_dir, ver, fcov)
