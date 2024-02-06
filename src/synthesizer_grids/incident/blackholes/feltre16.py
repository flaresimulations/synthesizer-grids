"""
Create a synthesizer incident grid for a broken power-law SED model using
the Feltre+2016 approach.
"""
import numpy as np
from unyt import c, Angstrom
from synthesizer_grids.parser import Parser
from synthesizer_grids.grid_io import GridFile


def broken_power_law(x, edges, indices, normalisations=False, normalise=True):
    """
    Function to produce a broken power-law spectrum.

    Arguments
        x (array, float)
            x-coordinate
        edges (list, float)
            List of floats describing the edges of the individual power-law
            components.
        indices (list, float)
            List of the power-law indices. Length should be one less than
            edges.
        normalisations (list, float)
            List of normalisations. Length should be one less than
            edges.
        normalise (bool)
            Flag to decide whether to normalise to unity.

    """

    y = np.zeros(len(x))

    # if normalisations not calculated, go ahead an
    if not normalisations:
        normalisations = [
            1.0,
        ]

        # calcualte remaining normalisations
        for i, (edge, ind1, ind2) in enumerate(
            zip(edges[1:], indices, indices[1:])
        ):
            norm = normalisations[i] * (edge**ind1) / (edge**ind2)
            normalisations.append(norm)

    # now construct spectra
    for e1, e2, ind, norm in zip(edges[0:],
                                 edges[1:],
                                 indices,
                                 normalisations):
        # identify indices within the wavelength range

        s = (x >= e1) & (x < e2)

        y[s] = norm * x[s] ** ind

    # normalise
    if normalise:
        y /= np.trapz(y[::-1], x=x[::-1])

    return y


if __name__ == "__main__":
    # Set up the command line arguments
    parser = Parser(description="Feltre et al. (2016) AGN model creation.")
    
    # path to grid directory (i.e. where incident and new grids are stored)
    # ** this can be replaced once Parser is updated. **
    parser.add_argument("-grid_dir", type=str, required=True)

    # get the arguments
    args = parser.parse_args()

    # Define the model name
    model = "feltre16"

    # Define the grid filename and path
    out_filename = f"{args.grid_dir}/{model}"

    alphas = np.arange(-2.0, -1.0, 0.2)
    axes = [alphas]

    # the shape of the grid (useful for creating outputs)
    axes_shape = list([len(axis) for axis in axes])

    # define edges
    edges_lam = [10.0, 2500.0, 100000.0, 1000000.0] * Angstrom  # Angstrom
    edges_nu = c / edges_lam
    edges = edges_nu.to("THz").value

    # define wavelength and frequency grid using the edges
    lam = np.arange(edges_lam.value[0], edges_lam.value[-1], 10.0) * Angstrom
    nu = c / lam
    x = nu.to("THz").value

    # define indices. The first entry is None because this is set by the axis.
    indices = [None, -0.5, 2.0]

    # create empty spectra grid
    spec = np.zeros((*axes_shape, len(lam)))

    # loop over alphas
    for i, alpha in enumerate(alphas):

        # use the value of alpha to fill in the missing power-law index
        indices[0] = alpha

        # generate the broken power-law spectra
        spec[i] = broken_power_law(x, edges[::-1], indices[::-1]) 


    # Create the GridFile ready to take outputs
    out_grid = GridFile(out_filename, mode="w", overwrite=True)

    # Write everything out thats common to all models
    out_grid.write_grid_common(
        model=model,
        axes={"alpha": alphas},
        wavelength=lam,
        spectra={"incident": spec * erg / s / Hz},
        alt_axes=("alphas"),
    )

    # Include the specific ionising photon luminosity
    print("Calculating and saving specific ionising luminosity")
    out_grid.add_specific_ionising_lum()