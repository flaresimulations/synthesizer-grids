"""
Create a synthesizer incident grid for the agnsed model
"""

import h5py
import numpy as np
from unyt import c, Angstrom

# from utils import __tag__, broken_power_law, add_specific_ionising_lum
from utils import broken_power_law, add_specific_ionising_lum
from datetime import date

# adding relagn to pythonpath
import sys

sys.path.append(
    "/Users/sw376/Dropbox/Research/projects/flares_emissiveagn/packages/RELAGN/src/"
)

from relagn import relagn

synthesizer_data_dir = "/Users/sw376/Dropbox/Research/data/synthesizer/"

model_name = "agnsed"


axes = ["mass", "accretion_rate_eddington", "cosine_inclination"]

axes_descriptions = {}
axes_descriptions["mass"] = "blackhole mass"
axes_descriptions[
    "accretion_rate_eddington"
] = "BH accretion rate / Edding accretion rate [LEdd=\eta MdotEdd c^2]"
axes_descriptions["cosine_inclination"] = "cosine of the inclination"

axes_units = {}
axes_units["mass"] = "Msun"
axes_units["accretion_rate_eddington"] = "dimensionless"
axes_units["cosine_inclination"] = "dimensionless"

axes_values = {}
axes_values["mass"] = 10 ** np.arange(6.0, 10.0, 1.0)
axes_values["accretion_rate_eddington"] = 10 ** np.arange(-3.0, 1.0, 1.0)
axes_values["cosine_inclination"] = [0.09, 0.5, 0.98]


# the shape of the grid (useful for creating outputs)
axes_shape = list([len(axes_values[axis]) for axis in axes])
print(axes_shape)

# filename
filename = f"{synthesizer_data_dir}/grids/dev/{model_name}.hdf5"

# open the new grid
with h5py.File(filename, "w") as hf:
    # save model attributes
    hf.attrs["model"] = model_name
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
    hf["spectra/wavelength"].attrs["Description"] = "Wavelength of the spectra grid"
    hf["spectra/wavelength"].attrs["Units"] = "Angstrom"

    # create empty spectra grid
    spec = np.zeros((*axes_shape, len(lam)))

    for i, M in enumerate(axes_values["mass"]):
        for j, accretion_rate_eddington in enumerate(
            axes_values["accretion_rate_eddington"]
        ):
            for k, cos_inc in enumerate(axes_values["cosine_inclination"]):
                dagn = relagn(
                    a=0.0,
                    cos_inc=cos_inc,
                    log_mdot=np.log10(accretion_rate_eddington),
                    M=M,
                )

                # lnu = dagn.get_totSED(rel=True) # relativistic
                lnu = dagn.get_totSED(rel=False)  # non-relativistic

                spec[i, j, k] = lnu[::-1]

    # save incident spectrum grid
    hf["spectra/incident"] = spec
    hf["spectra/incident"].attrs["Description"] = "Incident spectrum grid"
    hf["spectra/incident"].attrs["Units"] = "erg/s/Hz"


# calcualte specific_ionising_lum

add_specific_ionising_lum(filename)
