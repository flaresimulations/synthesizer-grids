import argparse
import numpy as np
import fsps
from datetime import date
from utils import (__tag__, write_data_h5py, write_attribute, add_log10Q, 
                   get_model_filename)

def generate_grid(model):
    """ Main function to create fsps grids used by synthesizer """

    

    # set up StellarPopulation grid
    imf_type = model['imf_type']
    imf_masses = model['imf_masses']

    if imf_type == 'chabrier03':
        sp = fsps.StellarPopulation(
            imf_type=1, imf_upper_limit=imf_masses[-1], imf_lower_limit=imf_masses[0])
    elif imf_type == 'bpl':
        imf_slopes = model['imf_slopes']
        # NOTE THIS ASSUMES boundaries of 0.5, 1.0

        if (imf_masses[1] != 0.5) or (imf_masses[2] != 1.0):
            # raise exception
            print(
                'WARNING: this IMF definition requires that the boundaries are [m_low, 0.5, 1.0, m_high]')

        sp = fsps.StellarPopulation(
            imf_type=2, imf_upper_limit=imf_masses[-1], imf_lower_limit=imf_masses[0], imf1=imf_slopes[0], imf2=imf_slopes[1], imf3=imf_slopes[2])
    else:
        # raise exception
        pass

    
    # include the isochrones/spectra as the sps_variant name
    # model['sps_variant'] = '-'.join(map(lambda x: str(x), sp.libraries[:2]))

    model['sps_variant'] = '-'.join(map(lambda x: x.decode("utf-8"), sp.libraries[:2]))
    print(model['sps_variant'])
    

    # generate the synthesizer_model_name
    synthesizer_model_name = get_model_filename(model)

    # this is the full path to the ultimate HDF5 grid file
    out_filename = f'{synthesizer_data_dir}/grids/{synthesizer_model_name}.hdf5'


    lam = sp.wavelengths  # units: Angstroms
    log10ages = sp.log_age  # units: log10(years)
    ages = 10**log10ages
    metallicities = sp.zlegend  # units: log10(Z)
    log10metallicities = np.log10(metallicities)

    na = len(log10ages)
    nZ = len(metallicities)

    log10Q = np.zeros((na, nZ))  # the ionising photon production rate
    spec = np.zeros((na, nZ, len(lam)))

    for iZ in range(nZ):
        spec_ = sp.get_spectrum(zmet=iZ+1)[1]   # 2D array Lsol / AA
        for ia in range(na):

            fnu = spec_[ia]  # Lsol / Hz
            fnu *= 3.826e33  # erg s^-1 Hz^-1 Msol^-1

            # fnu = convert_flam_to_fnu(lam, flam)
            spec[ia, iZ] = fnu

    # write out model parameters as top level attribute
    for key, value in model.items():
        write_attribute(out_filename, '/', key, (value))

    # write wavelength grid
    write_data_h5py(out_filename, 'spectra/wavelength', data=lam, overwrite=True)
    write_attribute(out_filename, 'spectra/wavelength', 'Description',
                    'Wavelength of the spectra grid')
    write_attribute(out_filename, 'spectra/wavelength', 'Units', 'Angstrom')

    # write stellar spectra
    write_data_h5py(out_filename, 'spectra/incident', data=spec, overwrite=True)
    write_attribute(out_filename, 'spectra/incident', 'Description',
                    'Three-dimensional spectra grid, [age, metallicity, wavelength]')
    write_attribute(out_filename, 'spectra/incident', 'Units', 'erg/s/Hz')

    # write out axes
    write_attribute(out_filename, '/', 'axes', ('log10age', 'metallicity'))

    # write out log10ages
    write_data_h5py(out_filename, 'axes/log10age', data=log10ages,
                    overwrite=True)
    write_attribute(out_filename, 'axes/log10age', 'Description',
                    'Stellar population ages in log10 years')
    write_attribute(out_filename, 'axes/log10age', 'Units', 'dex(yr)')

    # write out metallicities
    write_data_h5py(out_filename, 'axes/metallicity', data=metallicities, overwrite=True)
    write_attribute(out_filename, 'axes/metallicity', 'Description',
                    'raw abundances')
    write_attribute(out_filename, 'axes/metallicity', 'Units', 'dimensionless')

    return out_filename


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="FSPS download and grid creation")
    parser.add_argument('-synthesizer_data_dir', '--synthesizer_data_dir', default='/Users/sw376/Dropbox/Research/data/synthesizer')
    

    args = parser.parse_args()

    synthesizer_data_dir = args.synthesizer_data_dir
   
    grid_dir = f'{synthesizer_data_dir}/grids'
    
    default_model = {'sps_name': 'fsps',
                     
                     'sps_variant': False, # this is set later depending on the isochrones/spectra used
                     'imf_type': 'bpl',  # named IMF or bpl (broken power law)
                     'imf_masses': [0.08, 0.5, 1, 120],
                     'imf_slopes': [1.3, 2.3, 2.3],
                     'alpha': False,
                     'synthesizer-grids_tag': __tag__,
                     'date': str(date.today()), 
                     'pyfsps_version': str(fsps.__version__),
                     'sps_version': '3.2',
                     }

    models = []

    models += [{}]  # default model

    # # chabrier
    models += [{'imf_type': 'chabrier03', 'imf_masses': [0.08, 120], 'imf_slopes': []},  # chabrier03
               ]

    # different high-mass slopes
    models += [{'imf_slopes': [1.3, 2.3, a3]} for a3 in np.arange(2.8, 3.01, 0.1)]

    # different high-mass cut-offs
    models += [{'imf_type': 'chabrier03', 'imf_masses': [0.08, hmc]}
               for hmc in [1, 2, 5, 10, 20, 50, 100]]

    # different low-mass cut-offs
    models += [{'imf_type': 'chabrier03', 'imf_masses': [lmc, 120]}
               for lmc in [0.5, 1, 2, 5, 10, 20, 50]]

    for model_ in models:

        model = default_model | model_

        # make grid
        out_filename = generate_grid(model)

        # add log10Q for different ions
        add_log10Q(out_filename)
