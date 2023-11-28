"""
Generic incident grid creation functions module. 
"""

import os
import sys
import h5py
import numpy as np
from synthesizer.sed import calculate_Q
from synthesizer.photoionisation import Ions
from unyt import unyt_quantity

# import functions from grid_utils module
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from grid_utils import get_grid_properties_from_hdf5  # need this to work


def add_log10Q(grid_filename, ions=["HI", "HeII"], limit=100):
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
            hf[f"log10Q/{ion}"] = np.zeros(shape)

        # get spectra group
        spectra = hf[f"spectra"]

        # get wavelength grid, including with units
        # lam = spectra['wavelength'] * unyt_quantity.from_string(spectra['wavelength'].attrs['Units'])
        lam = spectra["wavelength"]

        # loop over grid points and calculate Q and store it
        for i, indices in enumerate(index_list):
            indices = tuple(indices)

            # loop over ions
            for ion in ions:
                # get the ionisation energy
                ionisation_energy = Ions.energy[ion]

                # get incident spectrum
                lnu = spectra["incident"][indices]

                # calculate Q
                Q = calculate_Q(
                    lam, lnu, ionisation_energy=ionisation_energy, limit=limit
                )

                # save
                hf[f"log10Q/{ion}"][indices] = np.log10(Q)
