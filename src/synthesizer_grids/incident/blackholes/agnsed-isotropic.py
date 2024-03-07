"""
Create a synthesizer incident grid for the agnsed model
"""

import h5py
import numpy as np
from unyt import c, Angstrom
from utils import (
    __tag__,
    broken_power_law,
    add_log10_specific_ionising_lum,
)
from datetime import date

# adding relagn to pythonpath
import sys

sys.path.append(
    "/Users/sw376/Dropbox/Research/projects/flares_emissiveagn/packages/RELAGN/src/"
)

from relagn import relagn

synthesizer_data_dir = "/Users/sw376/Dropbox/Research/data/synthesizer/"

model_name = "agnsed-isotropic"


axes = ["log10Mbh", "log10MdotEdd"]

axes_descriptions = {}
axes_descriptions["log10Mbh"] = "log10(BH mass)"
axes_descriptions[
    "log10MdotEdd"
] = "log10(BH accretion rate / Edding accretion rate) [LEdd=\eta MdotEdd c^2]"


axes_units = {}
axes_units["log10Mbh"] = "dex(Msun)"
axes_units["log10MdotEdd"] = "dimensionless"


axes_values = {}
axes_values["log10Mbh"] = np.arange(6.0, 10.0, 1.0)
axes_values["log10MdotEdd"] = np.arange(-3.0, 1.0, 1.0)


# the shape of the grid (useful for creating outputs)
axes_shape = list([len(axes_values[axis]) for axis in axes])
print(axes_shape)

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

    # initialise default model, to get wavelength grid
    dagn = relagn()
    lam = dagn.wave_grid[::-1]

    # save wavelength dataset
    hf["spectra/wavelength"] = lam
    hf["spectra/wavelength"].attrs[
        "Description"
    ] = "Wavelength of the spectra grid"
    hf["spectra/wavelength"].attrs["Units"] = "Angstrom"

    # create empty spectra grid
    spec = np.zeros((*axes_shape, len(lam)))

    for i, logM in enumerate(axes_values["log10Mbh"]):
        for j, log_mdot in enumerate(axes_values["log10MdotEdd"]):
            M = 10**logM

            dagn = relagn(a=0.0, log_mdot=log_mdot, M=M)

            # lnu = dagn.get_totSED(rel=True) # relativistic
            lnu = dagn.get_totSED(rel=False)  # non-relativistic

            spec[i, j] = lnu[::-1]

    # save incident spectrum grid
    hf["spectra/incident"] = spec
    hf["spectra/incident"].attrs["Description"] = "Incident spectrum grid"
    hf["spectra/incident"].attrs["Units"] = "erg/s/Hz"


# calcualte log10_specific_ionising_lum

add_log10_specific_ionising_lum(filename)
