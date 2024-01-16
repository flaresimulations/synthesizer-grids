"""
Download BPASS v2.3 and convert to HDF5 synthesizer grid.
"""
from hoki import load
import numpy as np
from synthesizer.sed import calc_log10_specific_ionising_lum
from synthesizer.cloudy import Ions
from datetime import date
from unyt import angstrom, erg, s, Hz

from synthesizer_grids.utilities.parser import Parser
from synthesizer_grids.utilities.grid_io import GridFile
from utils import get_model_filename


def resolve_name(original_model_name, bin, alpha=False):
    """Resolve the original BPASS model name into what we need. This is
    specific to 2.3. e.g. 'bpass_v2.3_chab300'"""

    bpass_imf = original_model_name.split("_")[-1]
    hmc = float(bpass_imf[-3:])  # high-mass cutoff

    if bpass_imf[:4] == "chab":
        imf_type = "chabrier03"
        imf_masses = [0.1, hmc]
        imf_slopes = False

    else:
        imf_type = "bpl"
        imf_masses = [0.1, 1.0, hmc]
        imf_slopes = [1.3, np.round(float(bpass_imf[:3]) / 100 + 1, 2)]

    model = {
        "original_model_name": original_model_name,
        "sps_name": "bpass",
        "sps_version": "2.3",
        "sps_variant": bin,
        "imf_type": imf_type,  # named IMF or bpl (broken power law)
        "imf_masses": imf_masses,
        "imf_slopes": imf_slopes,
        "alpha": alpha,
        "synthesizer-grids_tag": __tag__,
        "date": str(date.today()),
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


def make_single_alpha_grid(original_model_name, ae="+00", bs="bin"):
    """make a grid for a single alpha enhancement"""

    # convert bpass alpha code (e.g. '+02' into a numerical alpha e.g. 0.2)
    alpha = float(ae) / 10.0

    # returns a dictionary containing the sps model parameters
    model, bpass_imf = resolve_name(original_model_name, bs, alpha=alpha)

    # generate the synthesizer_model_name
    synthesizer_model_name = get_model_filename(model)

    print(synthesizer_model_name)

    # this is the full path to the ultimate HDF5 grid file
    out_filename = (
        f"{synthesizer_data_dir}/grids/{synthesizer_model_name}.hdf5"
    )

    # input directory
    input_dir = f'{synthesizer_data_dir}/input_files/bpass/{model["original_model_name"]}/'

    # create metallicity grid and dictionary
    map_key_to_met = {
        "zem5": 0.00001,
        "zem4": 0.0001,
        "z001": 0.001,
        "z002": 0.002,
        "z003": 0.003,
        "z004": 0.004,
        "z006": 0.006,
        "z008": 0.008,
        "z010": 0.01,
        "z014": 0.014,
        "z020": 0.020,
        "z030": 0.030,
        "z040": 0.040,
    }
    map_met_to_key = {k: v for v, k in map_key_to_met.items()}
    metallicities = np.sort(np.array(list(map_met_to_key.keys())))

    # get ages
    fn_ = f"{input_dir}/starmass-{bs}-imf_{bpass_imf}.a{ae}.{map_met_to_key[metallicities[0]]}.dat"
    starmass = load.model_output(fn_)
    log10ages = starmass["log_age"].values

    # get wavelength grid
    fn_ = f"spectra-{bs}-imf_{bpass_imf}.a{ae}.{map_met_to_key[metallicities[0]]}.dat"
    spec = load.model_output(f"{input_dir}/{fn_}")
    wavelengths = spec["WL"].values  # \AA
    nu = 3e8 / (wavelengths * 1e-10)

    # number of metallicities and ages
    nmetal = len(metallicities)
    na = len(log10ages)

    # set up outputs
    stellar_mass = np.zeros((na, nmetal))
    remnant_mass = np.zeros((na, nmetal))

    # the ionising photon production rate
    log10_specific_ionising_lum = {}
    log10_specific_ionising_lum["HI"] = np.zeros((na, nmetal))
    log10_specific_ionising_lum["HeII"] = np.zeros((na, nmetal))

    # provided by BPASS, sanity check for above
    log10_specific_ionising_lum_original = {}
    log10_specific_ionising_lum_original["HI"] = np.zeros((na, nmetal))

    spectra = np.zeros((na, nmetal, len(wavelengths)))

    for imetal, metal in enumerate(metallicities):
        print(imetal, metal)

        # get remaining and remnant fraction
        fn_ = f"{input_dir}/starmass-{bs}-imf_{bpass_imf}.a{ae}.{map_met_to_key[metal]}.dat"
        starmass = load.model_output(fn_)
        stellar_mass[:, imetal] = (
            starmass["stellar_mass"].values / 1e6
        )  # convert to per M_sol
        remnant_mass[:, imetal] = (
            starmass["remnant_mass"].values / 1e6
        )  # convert to per M_sol

        # get original log10_specific_ionising_lum
        fn_ = f"{input_dir}/ionizing-{bs}-imf_{bpass_imf}.a{ae}.{map_met_to_key[metal]}.dat"
        ionising = load.model_output(fn_)
        log10_specific_ionising_lum_original["HI"][:, imetal] = (
            ionising["prod_rate"].values - 6
        )  # convert to per M_sol

        # get spectra
        fn_ = f"{input_dir}/spectra-{bs}-imf_{bpass_imf}.a{ae}.{map_met_to_key[metal]}.dat"
        spec = load.model_output(fn_)

        for ia, log10age in enumerate(log10ages):
            spec_ = spec[str(log10age)].values  # Lsol AA^-1 10^6 Msol^-1

            # convert from Llam to Lnu
            spec_ /= 1e6  # Lsol AA^-1 Msol^-1
            spec_ *= 3.826e33  # erg s^-1 AA^-1 Msol^-1
            spec_ *= wavelengths / nu  # erg s^-1 Hz^-1 Msol^-1
            spectra[ia, imetal, :] = spec_

            # calcualte ionising photon luminosity
            for ion in ["HI", "HeII"]:
                limit = 100
                ionisation_energy = Ions.energy[ion]
                log10_specific_ionising_lum[ion][ia, imetal] = np.log10(
                    calc_log10_specific_ionising_lum(
                        wavelengths,
                        spec_,
                        ionisation_energy=ionisation_energy,
                        limit=limit,
                    )
                )
    # Create the GridFile ready to take outputs
    out_grid = GridFile(out_filename, mode="a", overwrite=True)

    # Write everything out thats common to all models
    out_grid.write_grid_common(
        model=model,
        axes={"log10age": log10ages, "metallicity": metallicities},
        wavelength=wavelengths * angstrom,
        spectra={"incident": spectra * erg / s / Hz},
        alt_axes=("log10ages", "metallicities"),
    )

    # Write datasets specific to BPASS 2.3
    out_grid.write_dataset(
        "star_fraction",
        stellar_mass,
        "Two-dimensional remaining stellar fraction grid, [age, Z]",
        units="Msun",
    )
    out_grid.write_dataset(
        "remnant_fraction",
        remnant_mass,
        "Two-dimensional remaining remnant fraction grid, [age, Z]",
        units="Msun",
    )
    out_grid.write_dataset(
        f"log10_specific_ionising_lum_original/{ion}",
        log10_specific_ionising_lum_original["HI"],
        "Two-dimensional (original) HI ionising photon production "
        "rate grid, [age,metal]",
        units="dimensionless",
    )

    # Include the specific ionising photon luminosity
    out_grid.add_specific_ionising_lum()


def make_full_grid(original_model_name, bs="bin"):
    """make a full grid for different alpha-ehancements"""

    # returns a dictionary containing the sps model parameters
    model, bpass_imf = resolve_name(original_model_name, bs)

    # generate the synthesizer_model_name
    synthesizer_model_name = get_model_filename(model)

    print(synthesizer_model_name)

    # this is the full path to the ultimate HDF5 grid file
    out_filename = (
        f"{synthesizer_data_dir}/grids/dev/{synthesizer_model_name}.hdf5"
    )

    # input directory
    input_dir = f'{synthesizer_data_dir}/input_files/bpass/{model["original_model_name"]}/'

    # --- ccreate metallicity grid and dictionary
    map_key_to_met = {
        "zem5": 0.00001,
        "zem4": 0.0001,
        "z001": 0.001,
        "z002": 0.002,
        "z003": 0.003,
        "z004": 0.004,
        "z006": 0.006,
        "z008": 0.008,
        "z010": 0.01,
        "z014": 0.014,
        "z020": 0.020,
        "z030": 0.030,
        "z040": 0.040,
    }
    map_met_to_key = {k: v for v, k in map_key_to_met.items()}
    metallicities = np.sort(np.array(list(map_met_to_key.keys())))
    log10metallicities = np.log10(metallicities)

    # --- create alpha-enhancement grid
    alpha_enhancements = np.array(
        [-0.2, 0.0, 0.2, 0.4, 0.6]
    )  # list of alpha enhancements
    ae_to_aek = {
        -0.2: "-02",
        0.0: "+00",
        0.2: "+02",
        0.4: "+04",
        0.6: "+06",
    }  # look up dictionary for filename

    # --- get ages
    fn_ = f"{input_dir}/starmass-bin-imf_{bpass_imf}.a+00.{map_met_to_key[metallicities[0]]}.dat"
    starmass = load.model_output(fn_)
    log10ages = starmass["log_age"].values

    # --- get wavelength grid
    fn_ = f"spectra-bin-imf_{bpass_imf}.a+00.{map_met_to_key[metallicities[0]]}.dat"
    spec = load.model_output(f"{input_dir}/{fn_}")
    wavelengths = spec["WL"].values  # \AA
    nu = 3e8 / (wavelengths * 1e-10)

    na = len(log10ages)
    nmetal = len(log10metallicities)
    nae = len(alpha_enhancements)

    # set up outputs
    stellar_mass = np.zeros((na, nmetal, nae))
    remnant_mass = np.zeros((na, nmetal, nae))

    # the ionising photon production rate
    log10_specific_ionising_lum = {}
    log10_specific_ionising_lum["HI"] = np.zeros((na, nmetal, nae))
    log10_specific_ionising_lum["HeII"] = np.zeros((na, nmetal, nae))

    # provided by BPASS, sanity check for above
    log10_specific_ionising_lum_original = {}
    log10_specific_ionising_lum_original["HI"] = np.zeros((na, nmetal, nae))

    spectra = np.zeros((na, nmetal, nae, len(wavelengths)))

    for imetal, metal in enumerate(metallicities):
        for iae, alpha_enhancement in enumerate(alpha_enhancements):
            print(metal, alpha_enhancement)

            aek = ae_to_aek[alpha_enhancement]
            metalk = map_met_to_key[metal]

            # --- get remaining and remnant fraction
            fn_ = f"{input_dir}/starmass-{bs}-imf_{bpass_imf}.a{aek}.{metalk}.dat"
            starmass = load.model_output(fn_)
            stellar_mass[:, imetal, iae] = (
                starmass["stellar_mass"].values / 1e6
            )
            remnant_mass[:, imetal, iae] = (
                starmass["remnant_mass"].values / 1e6
            )

            # --- get original log10_specific_ionising_lum
            fn_ = f"{input_dir}/ionizing-{bs}-imf_{bpass_imf}.a{aek}.{metalk}.dat"
            ionising = load.model_output(fn_)
            log10_specific_ionising_lum_original["HI"][:, imetal, iae] = (
                ionising["prod_rate"].values - 6
            )  # convert to per M_sol

            # --- get spectra
            fn_ = (
                f"{input_dir}/spectra-{bs}-imf_{bpass_imf}.a{aek}.{metalk}.dat"
            )
            spec = load.model_output(fn_)

            for ia, log10age in enumerate(log10ages):
                spec_ = spec[str(log10age)].values  # Lsol AA^-1 10^6 Msol^-1

                # --- convert from Llam to Lnu
                spec_ /= 1e6  # Lsol AA^-1 Msol^-1
                spec_ *= 3.826e33  # erg s^-1 AA^-1 Msol^-1
                spec_ *= wavelengths / nu  # erg s^-1 Hz^-1 Msol^-1

                spectra[ia, imetal, iae, :] = spec_  # Lsol AA^-1 10^6 Msol^-1

                # calcualte ionising photon luminosity
                for ion in ["HI", "HeII"]:
                    limit = 100
                    ionisation_energy = Ions.energy[ion]
                    log10_specific_ionising_lum[ion][
                        ia, imetal, iae
                    ] = np.log10(
                        calc_log10_specific_ionising_lum(
                            wavelengths,
                            spec_,
                            ionisation_energy=ionisation_energy,
                            limit=limit,
                        )
                    )

    # Create the GridFile ready to take outputs
    out_grid = GridFile(out_filename, mode="a", overwrite=True)

    # Write everything out thats common to all models
    out_grid.write_grid_common(
        model=model,
        axes={
            "log10age": log10ages,
            "metallicity": metallicities,
            "log10alpha": alpha_enhancements,
        },
        wavelength=wavelengths * angstrom,
        spectra={"incident": spectra * erg / s / Hz},
        alt_axes=("log10ages", "metallicities", "log10alphas"),
    )

    # Write datasets specific to BPASS 2.3
    out_grid.write_dataset(
        "star_fraction",
        stellar_mass,
        "Two-dimensional remaining stellar fraction grid, [age, Z]",
        units="Msun",
    )
    out_grid.write_dataset(
        "remnant_fraction",
        remnant_mass,
        "Two-dimensional remaining remnant fraction grid, [age, Z]",
        units="Msun",
    )
    out_grid.write_dataset(
        f"log10_specific_ionising_lum_original/{ion}",
        log10_specific_ionising_lum_original["HI"],
        "Two-dimensional (original) HI ionising photon production "
        "rate grid, [age,metal]",
        units="dimensionless",
    )

    # Include the specific ionising photon luminosity
    out_grid.add_specific_ionising_lum()


if __name__ == "__main__":
    # Set up the command line arguments
    parser = Parser(
        description="BPASS_2.3 download and grid creation",
        with_alpha=True,
    )
    parser.add_argument(
        "-models",
        "--models",
        default="bpass_v2.3_chab300",
        type=lambda arg: arg.split(","),
    )
    args = parser.parse_args()

    # Unpack the arguments
    synthesizer_data_dir = args.synthesizer_data_dir
    grid_dir = f"{synthesizer_data_dir}/grids"
    models = args.models

    print(models)

    for model in models:
        # The download currently doesn't work since these is no mirror
        # if args.download:
        #     download_data(model)
        #     untar_data(model)

        for bs in ["bin"]:  # no single star models , 'sin'
            # make a grid with a single alpha enahancement value
            if args.individual:
                for ae in ["-02", "+00", "+02", "+04", "+06"]:
                    # for ae in ['+00']: # used for testing
                    out_filename = make_single_alpha_grid(model, ae=ae, bs=bs)

            # make a full 3D grid
            if args.full:
                out_filename = make_full_grid(model, bs=bs)
