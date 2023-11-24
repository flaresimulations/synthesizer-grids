"""
Download Maraston2011 and convert to HDF5 synthesizer grid.
"""

#Allow the file to use incident_utils 

import numpy as np
import os
import argparse
from pathlib import Path
import tarfile
from synthesizer.utils import flam_to_fnu
from synthesizer.sed import calculate_Q
from datetime import date
import wget
import utils
import sys

#Allow the file to use incident_utils 
sys.path.append(os.path.join(os.path.dirname(__file__), '..')) 
from incident_utils import write_data_h5py, write_attribute, add_log10Q#, __tag__

c = 3*10**10


def download_data(input_dir):
    filename = wget.download(original_data_url)  # download the original data to the working directory

    Path(input_dir).mkdir(parents=True, exist_ok=True)

    # --- untar main directory
    tar = tarfile.open(filename)
    tar.extractall(path=input_dir)
    tar.close()
    os.remove(filename)



def make_grid(model, imf, extension):
    """Main function to convert Maraston 2011 and
    produce grids used by synthesizer"""

    #synthesizer_model_name = get_model_filename(model)
    #print(synthesizer_model_name)

    # Define output
    fname = f"{synthesizer_data_dir}/input_files/{model_name}/{model_name}{extension}_{imf}.hdf5" #ideally find a way to create the grids/model_name folder if it doesn't already exist

    metallicities = np.array([0.02])  # array of available metallicities

    log10metallicities = np.log10(metallicities)

    metallicity_code = {0.02: "002"}  # codes for converting metallicty

    fn = f"{input_dir}/ssp_M11_Pickles{extension}.{imf_code[imf]}z{metallicity_code[metallicities[0]]}"

    ages_, _, lam_, flam_ = np.loadtxt(fn).T # flam is in (ergs /s /AA /Msun)

    ages_Gyr = np.sort(np.array(list(set(ages_)))) # Gyr 
    ages = ages_Gyr * 1e9  # yr
    log10ages = np.log10(ages) 

    lam = lam_[ages_ == ages_[0]]

    flam_ = flam_ * (10**(-8))**2 / c #to convert from (ergs /s /AA /Msun) to (ergs /s /Hz /Msun)

    spec = np.zeros((len(ages), len(metallicities), len(lam)))

    for iZ, metallicity in enumerate(metallicities):
        for ia, age_Gyr in enumerate(ages_Gyr):
            print(iZ, ia, fn)
            ages_, _, lam_, flam_ = np.loadtxt(fn).T

            flam = flam_[ages_ == age_Gyr]
            fnu = flam_to_fnu(lam, flam)
            spec[ia, iZ] = fnu

    #write out spectra 
    write_data_h5py(fname, "spectra/wavelength", data=lam, overwrite=True)
    write_attribute(
        fname, "spectra/wavelength", "Description", "Wavelength of the spectra grid"
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


    #write out axes
    write_attribute(fname, '/', 'axes', ('log10age','metallicity'))

    write_data_h5py(fname, "axes/log10age", data=log10ages, overwrite=True)
    write_attribute(
        fname, "axes/log10age", "Description", "Stellar population ages in log10 years"
    )
    write_attribute(fname, "axes/log10age", "Units", "log10(yr)")

    write_data_h5py(fname, "axes/metallicity", data=metallicities, overwrite=True)
    write_attribute(fname, "axes/metallicity", "Description", "raw abundances")
    write_attribute(fname, "axes/metallicity", "Units", "dimensionless [Z]")

    return fname



# Lets include a way to call this script not via an entry point
if __name__ == "__main__":
    
    synthesizer_data_dir = '/home/sophienewman/code/synthesizer/src/synthesizer/data/synthesizer_data' #find smart way to do this

    model_name = "maraston11_pickles"

    input_dir = f"{synthesizer_data_dir}/original_data/{model_name}"  # the location to untar the original data

    original_data_url = "http://www.icg.port.ac.uk/~maraston/M11/SSP_M11_Pickles.tar.gz"
    imf_code = {"salpeter": "ss","kroupa": "kr","chabrier": "cha"}

    model = {'sps_name': 'maraston',
            'sps_version': False,
            'alpha': False,
            'date': str(date.today()), 
            } #'synthesizer-grids_tag': __tag__,

    for extension in ["","_nearIRextended","_UVtheoretical","_UVtheoretical_nearIRextended"]:
        
        if (extension == "") or (extension == "_nearIRextended"):
            imfs = ["salpeter","kroupa","chabrier"]

        #if imf == 'ss':
        #   model['imf_type'] = 'bpl'
        #   model['imf_masses'] = [0.1, 100]
        #   model['imf_slopes'] = [2.35]
        # this will work out the file name 

        if (extension == "_UVtheoretical") or (extension == "_UVtheoretical_nearIRextended"):
            imfs = ["salpeter"]

        for imf in imfs: 

            print(extension)
            download_data(input_dir)
            fname = make_grid(model,imf,extension) #makes the grid and returns the name 

            add_log10Q(fname)

            

    