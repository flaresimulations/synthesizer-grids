"""
Create a synthesizer incident grid for a broken power-law SED model.
"""

import h5py 
import numpy as np
from unyt import c, Angstrom
from utils import __tag__, broken_power_law, add_log10Q
from datetime import date

synthesizer_data_dir = '/Users/sw376/Dropbox/Research/data/synthesizer/'

model_name = 'feltre16'


axes = ['alpha']

axes_descriptions = {}
axes_descriptions['alpha'] = 'X-ray - UV power-law slope'

axes_units = {}
axes_units['alpha'] = 'dimensionless'

axes_values = {}
axes_values['alpha'] = np.arange(-2.0, -1.0, 0.2)


# the shape of the grid (useful for creating outputs)
axes_shape = list([len(axes_values[axis]) for axis in axes])
print(axes_shape)

# define edges
edges_lam = [10.0, 2500.0, 100000.0, 1000000.0] * Angstrom  # Angstrom
edges_nu = c / edges_lam
edges = edges_nu.to('Hz').value

# define wavelength and frequency grid
lam = np.arange(edges_lam.value[0], edges_lam.value[-1], 10.) * Angstrom
nu = c / lam
x = nu.to('Hz').value

# define indices
indices = [None, -0.5, 2.0]


# filename
filename = f'{synthesizer_data_dir}/grids/{model_name}.hdf5'

# open the new grid
with h5py.File(filename, 'w') as hf:

    # save model attributes
    hf.attrs['model'] = model_name
    hf.attrs['synthesizer-grids_tag'] = __tag__
    hf.attrs['date'] = str(date.today())
    hf.attrs['axes'] = axes

    # save axis quantities
    for axis in axes:
        hf[f'axes/{axis}'] = axes_values[axis]
        hf[f'axes/{axis}'].attrs['Description'] = axes_descriptions[axis]
        hf[f'axes/{axis}'].attrs['Units'] = axes_units[axis]


    # create empty spectra grid
    spec = np.zeros((*axes_shape, len(lam)))
    
    for i, alpha in enumerate(axes_values['alpha']):

        indices[0] = alpha
        spec[i] = broken_power_law(x, edges, indices)


    # save wavelength dataset
    hf['spectra/wavelength'] = lam
    hf['spectra/wavelength'].attrs['Description'] = 'Wavelength of the spectra grid'
    hf['spectra/wavelength'].attrs['Units'] = 'Angstrom'

    # save incident spectrum grid
    hf['spectra/incident'] = spec
    hf['spectra/incident'].attrs['Description'] = 'Incident spectrum grid'
    hf['spectra/incident'].attrs['Units'] = 'erg/s/Hz'


# calcualte log10Q

add_log10Q(filename)
