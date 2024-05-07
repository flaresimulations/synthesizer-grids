import sys
import os
import numpy as np
from unyt import c
import scipy.integrate as integrate

# import function from incident_utils module
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
import incident_utils
from incident_utils import *

# __tag__ = incident_utils.__tag__


def broken_power_law(x, edges, indices, norms=False, normalise=True):
    """
    Function to produce a broken power-law spectrum

    """

    y = np.zeros(len(x))

    # if normalisations not calculated, go ahead an
    if not norms:
        norms = [
            1.0,
        ]

        # calcualte remaining normalisations
        for i, (edge, ind1, ind2) in enumerate(
            zip(edges[1:], indices, indices[1:])
        ):
            norm = norms[i] * (edge**ind1) / (edge**ind2)
            norms.append(norm)

    # now construct spectra
    for e1, e2, ind, norm in zip(edges[0:], edges[1:], indices, norms):
        # identify indices within the wavelength range

        s = (x >= e1) & (x < e2)

        print(e1, e2, ind, norm, np.sum(s))
        y[s] = norm * x[s] ** ind

    # normalise
    if normalise:
        y /= np.trapz(y[::-1], x=x[::-1])

    return y
