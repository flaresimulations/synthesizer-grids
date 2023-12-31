"""
Download BPASS v2.3 and convert to HDF5 synthesizer grid.
"""

from hoki import load
import argparse
import numpy as np
import gdown
import tarfile
import h5py
from scipy import integrate
from unyt import h, c
from synthesizer.sed import calculate_Q
from synthesizer.cloudy import Ions
from datetime import date
from utils import (__tag__, write_data_h5py, write_attribute, add_log10Q, 
                   get_model_filename)

def resolve_name(original_model_name, bin, alpha = False):
    """ Resolve the original BPASS model name into what we need. This is specific to 2.3. e.g. 'bpass_v2.3_chab300' """

    bpass_imf = original_model_name.split('_')[-1]
    hmc = float(bpass_imf[-3:])  # high-mass cutoff

    if bpass_imf[:4] == 'chab':
        imf_type = 'chabrier03'
        imf_masses = [0.1, hmc]
        imf_slopes = False
    
    else:
        imf_type = 'bpl'
        imf_masses = [0.1, 1.0, hmc]
        imf_slopes = [1.3, np.round(float(bpass_imf[:3])/100+1, 2)]

    model = {'original_model_name': original_model_name,
             'sps_name': 'bpass',
             'sps_version': '2.3',
             'sps_variant': bin,
             'imf_type': imf_type,  # named IMF or bpl (broken power law)
             'imf_masses': imf_masses,
             'imf_slopes': imf_slopes,
             'alpha': alpha,
             'synthesizer-grids_tag': __tag__,
             'date': str(date.today()), 
             }

    print(model)

    return model, bpass_imf


# # not currently used
# def download_data(model):

#     if model in model_url.keys():
#         filename = gdown.download(model_url[model], quiet=False, fuzzy=True)
#         return filename
#     else:
#         print('ERROR: no url for that model')

# # not currently used
# def untar_data(model, remove_archive = False):


#     input_dir = f'{parent_model_dir}/{model}'
#     tar = tarfile.open(f'{parent_model_dir}/{model}.tar')
#     tar.extractall(path = input_dir)
#     tar.close()
#     if remove_archive: os.remove(f'{parent_model_dir}/{model}.tar')


def make_single_alpha_grid(original_model_name, ae = '+00', bs = 'bin'):

    """ make a grid for a single alpha enhancement """

    # convert bpass alpha code (e.g. '+02' into a numerical alpha e.g. 0.2)
    alpha = float(ae)/10.

    # returns a dictionary containing the sps model parameters
    model, bpass_imf = resolve_name(original_model_name, bs, alpha = alpha)

    # generate the synthesizer_model_name
    synthesizer_model_name = get_model_filename(model)

    print(synthesizer_model_name)

    # this is the full path to the ultimate HDF5 grid file
    out_filename = f'{synthesizer_data_dir}/grids/{synthesizer_model_name}.hdf5'

    # input directory
    input_dir = f'{synthesizer_data_dir}/input_files/bpass/{model["original_model_name"]}/'

    # create metallicity grid and dictionary
    Zk_to_Z = {'zem5': 0.00001, 'zem4': 0.0001, 'z001': 0.001, 'z002': 0.002, 'z003': 0.003, 'z004': 0.004, 'z006': 0.006, 'z008': 0.008, 'z010': 0.01, 'z014': 0.014, 'z020': 0.020, 'z030': 0.030, 'z040': 0.040}
    Z_to_Zk = {k:v for v, k in Zk_to_Z.items()}
    Zs = np.sort(np.array(list(Z_to_Zk.keys())))
    
    # get ages
    fn_ = f'{input_dir}/starmass-{bs}-imf_{bpass_imf}.a{ae}.{Z_to_Zk[Zs[0]]}.dat'
    starmass = load.model_output(fn_)
    log10ages = starmass['log_age'].values

    # get wavelength grid
    fn_ = f'spectra-{bs}-imf_{bpass_imf}.a{ae}.{Z_to_Zk[Zs[0]]}.dat'
    spec = load.model_output(f'{input_dir}/{fn_}')
    wavelengths = spec['WL'].values # \AA
    nu = 3E8/(wavelengths*1E-10)

    # number of metallicities and ages
    nZ = len(Zs)
    na = len(log10ages)

    # set up outputs
    stellar_mass = np.zeros((na,nZ))
    remnant_mass = np.zeros((na,nZ))

    # the ionising photon production rate
    log10Q = {}
    log10Q['HI'] = np.zeros((na,nZ)) 
    log10Q['HeII'] = np.zeros((na,nZ)) 

    # provided by BPASS, sanity check for above
    log10Q_original = {}
    log10Q_original['HI'] = np.zeros((na,nZ)) 
    
    spectra = np.zeros((na, nZ, len(wavelengths)))

    for iZ, Z in enumerate(Zs):

        print(iZ, Z)

        # get remaining and remnant fraction
        fn_ = f'{input_dir}/starmass-{bs}-imf_{bpass_imf}.a{ae}.{Z_to_Zk[Z]}.dat'
        starmass = load.model_output(fn_)
        stellar_mass[:, iZ] = starmass['stellar_mass'].values/1E6  # convert to per M_sol
        remnant_mass[:, iZ] = starmass['remnant_mass'].values/1E6  # convert to per M_sol

        # get original log10Q
        fn_ = f'{input_dir}/ionizing-{bs}-imf_{bpass_imf}.a{ae}.{Z_to_Zk[Z]}.dat'
        ionising = load.model_output(fn_)
        log10Q_original['HI'][:, iZ] = ionising['prod_rate'].values - 6   # convert to per M_sol

        # get spectra
        fn_ = f'{input_dir}/spectra-{bs}-imf_{bpass_imf}.a{ae}.{Z_to_Zk[Z]}.dat'
        spec = load.model_output(fn_)

        for ia, log10age in enumerate(log10ages):

            spec_ = spec[str(log10age)].values # Lsol AA^-1 10^6 Msol^-1

            # convert from Llam to Lnu
            spec_ /= 1E6 # Lsol AA^-1 Msol^-1
            spec_ *= (3.826e33)  # erg s^-1 AA^-1 Msol^-1
            spec_ *= wavelengths/nu # erg s^-1 Hz^-1 Msol^-1
            spectra[ia, iZ, :] = spec_

            # calcualte ionising photon luminosity
            for ion in ['HI', 'HeII']:
                limit = 100
                ionisation_energy = Ions.energy[ion]
                log10Q[ion][ia, iZ] = np.log10(calculate_Q(wavelengths, spec_, ionisation_energy=ionisation_energy, limit=limit))


    # write out model parameters as top level attribute
    for key, value in model.items():
        # print(key, value)
        write_attribute(out_filename, '/', key, (value))

    write_data_h5py(out_filename, 'star_fraction', data=stellar_mass, overwrite=True)
    write_attribute(out_filename, 'star_fraction', 'Description',
                    'Two-dimensional remaining stellar fraction grid, [age,Z]')

    write_data_h5py(out_filename, 'remnant_fraction', data=remnant_mass, overwrite=True)
    write_attribute(out_filename, 'remnant_fraction', 'Description',
                    'Two-dimensional remaining remnant fraction grid, [age,Z]')

    for ion in ['HI']:
        write_data_h5py(out_filename, f'log10Q_original/{ion}', data=log10Q_original[ion], overwrite=True)
        write_attribute(out_filename, f'log10Q_original/{ion}', 'Description',
                        f'Two-dimensional (original) {ion} ionising photon production rate grid, [age,Z]')
        write_attribute(out_filename, f'log10Q_original/{ion}', 'Units',
                        'dex(1/s)')

    for ion in ['HI', 'HeII']:
        write_data_h5py(out_filename, f'log10Q/{ion}', data=log10Q[ion], overwrite=True)
        write_attribute(out_filename, f'log10Q/{ion}', 'Description',
                        f'Two-dimensional {ion} ionising photon production rate grid, [age,Z]')
        write_attribute(out_filename, f'log10Q/{ion}', 'Units',
                        'dex(1/s)')

    write_data_h5py(out_filename, 'spectra/incident', data=spectra, overwrite=True)
    write_attribute(out_filename, 'spectra/incident', 'Description',
                    'Three-dimensional spectra grid, [Z,Age,wavelength]')
    write_attribute(out_filename, 'spectra/incident', 'Units', 'erg/s/Hz')

    write_data_h5py(out_filename, 'spectra/wavelength', data=wavelengths, overwrite=True)
    write_attribute(out_filename, 'spectra/wavelength', 'Description',
            'Wavelength of the spectra grid')
    write_attribute(out_filename, 'spectra/wavelength', 'Units', 'Angstrom')


    # write out axes
    write_attribute(out_filename, '/', 'axes', ('log10age', 'metallicity'))

    # write out log10ages
    write_data_h5py(out_filename, 'axes/log10age', data=log10ages,
                    overwrite=True)
    write_attribute(out_filename, 'axes/log10age', 'Description',
                    'Stellar population ages in log10 years')
    write_attribute(out_filename, 'axes/log10age', 'Units', 'dex(yr)')

    # write out metallicities
    write_data_h5py(out_filename, 'axes/metallicity', data=Zs, overwrite=True)
    write_attribute(out_filename, 'axes/metallicity', 'Description',
                    'raw abundances')
    write_attribute(out_filename, 'axes/metallicity', 'Units', 'dimensionless')

    return out_filename


def make_full_grid(original_model_name, bs = 'bin'):

    """ make a full grid for different alpha-ehancements """

    # returns a dictionary containing the sps model parameters
    model, bpass_imf = resolve_name(original_model_name, bs)

    # generate the synthesizer_model_name
    synthesizer_model_name = get_model_filename(model)

    print(synthesizer_model_name)

    # this is the full path to the ultimate HDF5 grid file
    out_filename = f'{synthesizer_data_dir}/grids/dev/{synthesizer_model_name}.hdf5'

    # input directory
    input_dir = f'{synthesizer_data_dir}/input_files/bpass/{model["original_model_name"]}/'

    # --- ccreate metallicity grid and dictionary
    Zk_to_Z = {'zem5': 0.00001, 'zem4': 0.0001, 'z001': 0.001, 'z002': 0.002, 'z003': 0.003, 'z004': 0.004, 'z006': 0.006, 'z008': 0.008, 'z010': 0.01, 'z014': 0.014, 'z020': 0.020, 'z030': 0.030, 'z040': 0.040}
    Z_to_Zk = {k:v for v, k in Zk_to_Z.items()}
    Zs = np.sort(np.array(list(Z_to_Zk.keys())))
    log10Zs = np.log10(Zs)

    # --- create alpha-enhancement grid
    alpha_enhancements = np.array([-0.2, 0.0, 0.2, 0.4, 0.6])  # list of alpha enhancements
    ae_to_aek = {-0.2:'-02', 0.0:'+00', 0.2:'+02', 0.4:'+04', 0.6:'+06'}  # look up dictionary for filename

    # --- get ages
    fn_ = f'{input_dir}/starmass-bin-imf_{bpass_imf}.a+00.{Z_to_Zk[Zs[0]]}.dat'
    starmass = load.model_output(fn_)
    log10ages = starmass['log_age'].values

    # --- get wavelength grid
    fn_ = f'spectra-bin-imf_{bpass_imf}.a+00.{Z_to_Zk[Zs[0]]}.dat'
    spec = load.model_output(f'{input_dir}/{fn_}')
    wavelengths = spec['WL'].values # \AA
    nu = 3E8/(wavelengths*1E-10)

    na = len(log10ages)
    nZ = len(log10Zs)
    nae = len(alpha_enhancements)



    # set up outputs
    stellar_mass = np.zeros((na,nZ,nae))
    remnant_mass = np.zeros((na,nZ,nae))

    # the ionising photon production rate
    log10Q = {}
    log10Q['HI'] = np.zeros((na,nZ,nae)) 
    log10Q['HeII'] = np.zeros((na,nZ,nae)) 

    # provided by BPASS, sanity check for above
    log10Q_original = {}
    log10Q_original['HI'] = np.zeros((na,nZ,nae)) 


    spectra = np.zeros((na, nZ, nae, len(wavelengths)))



    for iZ, Z in enumerate(Zs):

        for iae, alpha_enhancement in enumerate(alpha_enhancements):

            print(Z, alpha_enhancement)

            aek = ae_to_aek[alpha_enhancement]
            Zk = Z_to_Zk[Z]

            # --- get remaining and remnant fraction
            fn_ = f'{input_dir}/starmass-{bs}-imf_{bpass_imf}.a{aek}.{Zk}.dat'
            starmass = load.model_output(fn_)
            stellar_mass[:, iZ, iae] = starmass['stellar_mass'].values/1E6
            remnant_mass[:, iZ, iae] = starmass['remnant_mass'].values/1E6

            # --- get original log10Q
            fn_ = f'{input_dir}/ionizing-{bs}-imf_{bpass_imf}.a{aek}.{Zk}.dat'
            ionising = load.model_output(fn_)
            log10Q_original['HI'][:, iZ, iae] = ionising['prod_rate'].values - 6   # convert to per M_sol

            # --- get spectra
            fn_ = f'{input_dir}/spectra-{bs}-imf_{bpass_imf}.a{aek}.{Zk}.dat'
            spec = load.model_output(fn_)

            for ia, log10age in enumerate(log10ages):

                spec_ = spec[str(log10age)].values # Lsol AA^-1 10^6 Msol^-1

                # --- convert from Llam to Lnu
                spec_ /= 1E6 # Lsol AA^-1 Msol^-1
                spec_ *= (3.826e33)  # erg s^-1 AA^-1 Msol^-1
                spec_ *= wavelengths/nu # erg s^-1 Hz^-1 Msol^-1

                spectra[ia, iZ, iae, :] = spec_ # Lsol AA^-1 10^6 Msol^-1

                # calcualte ionising photon luminosity
                for ion in ['HI', 'HeII']:
                    limit = 100
                    ionisation_energy = Ions.energy[ion]
                    log10Q[ion][ia, iZ, iae] = np.log10(calculate_Q(wavelengths, spec_, ionisation_energy=ionisation_energy, limit=limit))


    # write out model parameters as top level attribute
    for key, value in model.items():
        # print(key, value)
        write_attribute(out_filename, '/', key, (value))

    write_data_h5py(out_filename, 'star_fraction', data=stellar_mass, overwrite=True)
    write_attribute(out_filename, 'star_fraction', 'Description',
                    'Two-dimensional remaining stellar fraction grid, [age,Z]')

    write_data_h5py(out_filename, 'remnant_fraction', data=remnant_mass, overwrite=True)
    write_attribute(out_filename, 'remnant_fraction', 'Description',
                    'Two-dimensional remaining remnant fraction grid, [age,Z]')


    # write out ionising photon production rate
    for ion in ['HI']:
        write_data_h5py(out_filename, f'log10Q_original/{ion}', data=log10Q_original[ion], overwrite=True)
        write_attribute(out_filename, f'log10Q_original/{ion}', 'Description',
                        f'Two-dimensional (original) {ion} ionising photon production rate grid, [age,Z]')
        write_attribute(out_filename, f'log10Q_original/{ion}', 'Units',
                        'dex(1/s)')

    for ion in ['HI', 'HeII']:
        write_data_h5py(out_filename, f'log10Q/{ion}', data=log10Q[ion], overwrite=True)
        write_attribute(out_filename, f'log10Q/{ion}', 'Description',
                        f'Two-dimensional {ion} ionising photon production rate grid, [age,Z]')
        write_attribute(out_filename, f'log10Q/{ion}', 'Units',
                        'dex(1/s)')


    # write out axes
    write_attribute(out_filename, '/', 'axes', ('log10age', 'metallicity', 'alpha'))

    # write out log10ages
    write_data_h5py(out_filename, 'axes/log10age', data=log10ages,
                    overwrite=True)
    write_attribute(out_filename, 'axes/log10age', 'Description',
                    'Stellar population ages in log10 years')
    write_attribute(out_filename, 'axes/log10age', 'Units', 'dex(yr)')

    # write out metallicities
    write_data_h5py(out_filename, 'axes/metallicity', data=Zs, overwrite=True)
    write_attribute(out_filename, 'axes/metallicity', 'Description',
                    'raw abundances')
    write_attribute(out_filename, 'axes/metallicity', 'Units', 'dimensionless')

    # write of alpha values
    write_data_h5py(out_filename, 'axes/alpha', data=alpha_enhancements, overwrite=True)
    write_attribute(out_filename, 'axes/alpha', 'Description',
                    'log10(alpha enhancement)')
    write_attribute(out_filename, 'axes/alpha', 'Units', 'dimensionless')




    # write wavelength grid
    write_data_h5py(out_filename, 'spectra/wavelength', data=wavelengths, overwrite=True)
    write_attribute(out_filename, 'spectra/wavelength', 'Description',
            'Wavelength of the spectra grid')
    write_attribute(out_filename, 'spectra/wavelength', 'Units', 'Angstrom')

    # write stellar spectra 
    write_data_h5py(out_filename, 'spectra/stellar', data=spectra, overwrite=True)
    write_attribute(out_filename, 'spectra/stellar', 'Description',
                    'Three-dimensional spectra grid, [age, Z, ae, wavelength]')
    write_attribute(out_filename, 'spectra/stellar', 'Units', 'erg/s/Hz')




    return out_filename










if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="BPASS_2.3 download and grid creation")
    
    # flag whether to download data
    parser.add_argument('--download', default=False, action='store_true',
                        help=("download bpass data directly in current directory "
                              "and untar in sunthesizer data directory"))

    # path to synthesizer dir
    parser.add_argument('-synthesizer_data_dir', '--synthesizer_data_dir', default='/Users/sw376/Dropbox/Research/data/synthesizer')

    # flag whether to make individual alpha grids
    parser.add_argument('-individual', '--individual', action=argparse.BooleanOptionalAction)

    # flag whether to make full grid
    parser.add_argument('-full', '--full', action=argparse.BooleanOptionalAction, default=False)

    # models
    parser.add_argument('-models', '--models', default='bpass_v2.3_chab300')

    # arguments    
    args = parser.parse_args()

    
    synthesizer_data_dir = args.synthesizer_data_dir

    # get grid dir
    grid_dir = f'{synthesizer_data_dir}/grids'

    models = args.models.split(',')

    print(models)

    for model in models:

        # if args.download:
        #     download_data(model)
        #     untar_data(model)

        for bs in ['bin']: # no single star models , 'sin'

            # make a grid with a single alpha enahancement value
            if args.individual:

                for ae in ['-02','+00','+02','+04','+06']:
                # for ae in ['+00']: # used for testing
                    out_filename = make_single_alpha_grid(model, ae = ae, bs = bs)

            # make a full 3D grid 
            if args.full:

                out_filename = make_full_grid(model, bs = bs)
