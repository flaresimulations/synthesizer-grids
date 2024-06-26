"""
Generic incident grid creation functions module.
"""

import os
import sys

import h5py
import numpy as np
from synthesizer.photoionisation import Ions
from synthesizer.sed import Sed
from unyt import Angstrom

# import functions from grid_utils module
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from grid_utils import get_grid_properties_from_hdf5  # need this to work

# __tag__ = grid_utils.__tag__


def add_log10_specific_ionising_lum(
    grid_filename, ions=["HI", "HeII"], limit=100
):
    """
    A function to calculate the ionising photon luminosity for different ions.

    Parameters
    ---------
    grid_filename : str
        the filename of the HDF5 grid
    ions : list
        a list of ions to calculate Q for
    limit: float or int, optional
        An upper bound on the number of subintervals
        used in the integration adaptive algorithm.

    """
    # open  new grid
    with h5py.File(grid_filename, "a") as hf:
        # Get the properties of the grid including the dimensions etc.
        (
            axes,
            shape,
            n_models,
            mesh,
            model_list,
            index_list,
        ) = get_grid_properties_from_hdf5(hf)

        # set up output arrays
        for ion in ions:
            hf[f"log10_specific_ionising_lum/{ion}"] = np.zeros(shape)

        # get spectra group
        spectra = hf["spectra"]

        # get wavelength grid, including with units
        lam = np.array(spectra["wavelength"]) * Angstrom

        # loop over grid points and calculate Q and store it
        for i, indices in enumerate(index_list):
            indices = tuple(indices)

            # loop over ions
            for ion in ions:
                # get the ionisation energy
                ionisation_energy = Ions.energy[ion]

                # get incident spectrum
                lnu = spectra["incident"][indices]

                # create sed object
                sed = Sed(lam=lam, lnu=lnu)

                # calculate Q
                Q = sed.calculate_ionising_photon_production_rate(
                    ionisation_energy=ionisation_energy, limit=limit
                )

                # save
                hf[f"log10_specific_ionising_lum/{ion}"][indices] = np.log10(Q)


def write_data_h5py(filename, name, data, overwrite=False):
    check = check_h5py(filename, name)

    with h5py.File(filename, "a") as h5file:
        if check:
            if overwrite:
                print("Overwriting data in %s" % name)
                del h5file[name]
                h5file[name] = data
            else:
                raise ValueError(
                    "Dataset already exists, " + "and `overwrite` not set"
                )
        else:
            h5file.create_dataset(name, data=data)


def check_h5py(filename, obj_str):
    with h5py.File(filename, "a") as h5file:
        if obj_str not in h5file:
            return False
        else:
            return True


def load_h5py(filename, obj_str):
    with h5py.File(filename, "a") as h5file:
        dat = np.array(h5file.get(obj_str))
    return dat


def write_attribute(filename, obj, key, value):
    """
    Write attribute to an HDF5 file

    Args
    obj (str) group  or dataset to attach attribute to
    key (str) attribute key string
    value (str) content of the attribute
    """
    with h5py.File(filename, "a") as h5file:
        dset = h5file[obj]
        dset.attrs[key] = value


def get_names_h5py(filename, group):
    """
    Return list of the names of objects inside a group
    """
    with h5py.File(filename, "r") as h5file:
        keys = list(h5file[group].keys())

    return keys


def load_arr(name, filename):
    """
    Load Dataset array from file
    """
    with h5py.File(filename, "r") as f:
        if name not in f:
            raise ValueError("'%s' Dataset doesn't exist..." % name)

        arr = np.array(f.get(name))

    return arr


def read_params(param_file):
    """
    Args:
    param_file (str) location of parameter file

    Returns:
    parameters (object)
    """
    return __import__(param_file)
