"""
Create a synthesizer incident grid for a broken power-law SED model.
"""

from datetime import date

import h5py
import numpy as np
from unyt import Angstrom, c
from utils import (
    __tag__,
    add_log10_specific_ionising_lum,
    broken_power_law,
)

synthesizer_data_dir = "/Users/sw376/Dropbox/Research/data/synthesizer/"

model_name = "feltre16"


# axes = ['alpha', 'metallicity']
axes = ["alpha"]

axes_descriptions = {}
axes_descriptions["alpha"] = "X-ray - UV power-law slope"
# axes_descriptions['metallicity'] = 'gas phase metallicity'

axes_units = {}
axes_units["alpha"] = "dimensionless"
# axes_units['metallicity'] = 'dimensionless'

axes_values = {}
axes_values["alpha"] = np.arange(-2.0, -1.0, 0.2)
# axes_values['metallicity'] = [0.00001, 0.00003, 0.0001,
#                               0.0003, 0.001, 0.003, 0.01, 0.03]
# axes_values['alpha'] = [-2.0]

# the shape of the grid (useful for creating outputs)
axes_shape = list([len(axes_values[axis]) for axis in axes])
print(axes_shape)

# define edges
edges_lam = [10.0, 2500.0, 100000.0, 1000000.0] * Angstrom  # Angstrom
edges_nu = c / edges_lam
edges = edges_nu.to("THz").value

# define wavelength and frequency grid
lam = np.arange(edges_lam.value[0], edges_lam.value[-1], 10.0) * Angstrom
nu = c / lam
x = nu.to("THz").value

# define indices
indices = [None, -0.5, 2.0]


# filename
filename = f"{synthesizer_data_dir}/grids/dev/{model_name}.hdf5"

# open the new grid
with h5py.File(filename, "w") as hf:
    # save model attributes
    hf.attrs["model"] = model_name
    hf.attrs["synthesizer-grids_tag"] = __tag__
    hf.attrs["date"] = str(date.today())
    hf.attrs["axes"] = axes

    # save axis quantities
    for axis in axes:
        hf[f"axes/{axis}"] = axes_values[axis]
        hf[f"axes/{axis}"].attrs["Description"] = axes_descriptions[axis]
        hf[f"axes/{axis}"].attrs["Units"] = axes_units[axis]

    # create empty spectra grid
    spec = np.zeros((*axes_shape, len(lam)))

    for i, alpha in enumerate(axes_values["alpha"]):
        indices[0] = alpha
        lnu = broken_power_law(x, edges[::-1], indices[::-1])
        spec[i] = lnu

        # for j, metallicity in enumerate(axes_values['metallicity']):
        #     spec[i, j] = lnu

    # save wavelength dataset
    hf["spectra/wavelength"] = lam
    hf["spectra/wavelength"].attrs["Description"] = (
        "Wavelength of the spectra grid"
    )
    hf["spectra/wavelength"].attrs["Units"] = "Angstrom"

    # save incident spectrum grid
    hf["spectra/incident"] = spec
    hf["spectra/incident"].attrs["Description"] = "Incident spectrum grid"
    hf["spectra/incident"].attrs["Units"] = "erg/s/Hz"


# calcualte log10_specific_ionising_lum

add_log10_specific_ionising_lum(filename)
