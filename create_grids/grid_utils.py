
import numpy as np
# import git

#repo = git.Repo(search_parent_directories=True)
#__tag__ = repo.git.describe()
#__tag__ = __version__
#print(__tag__)



def get_grid_properties_from_hdf5(hf, verbose=True):

    """
    A wrapper over get_grid_properties to get the grid properties for a HDF5
    grid.
    """
    axes = hf.attrs['axes']  # list of axes

    # dictionary of axis grid points
    axes_values = {axis: hf['axes'][axis][:] for axis in axes}
    
    # Get the properties of the grid including the dimensions etc.
    return get_grid_properties(axes, axes_values, verbose=verbose)


def get_grid_properties(axes, axes_values, verbose=True):

    """ 
    Get the properties of the grid including the dimensions etc.
    """

    # the grid axes   
    if verbose: print(f'axes: {axes}')

    # number of axes
    n_axes = len(axes)
    if verbose: print(f'number of axes: {n_axes}')

    # the shape of the grid (useful for creating outputs)
    shape = list([len(axes_values[axis]) for axis in axes])
    if verbose: print(f'shape: {shape}')

    # determine number of models
    n_models = np.prod(shape)
    if verbose: print(f'number of models to run: {n_models}')

    # create the mesh of the grid
    mesh = np.array(np.meshgrid(*[np.array(axes_values[axis]) for axis in axes]))

    # create the list of the models 
    model_list = mesh.T.reshape(n_models, n_axes)
    if verbose: 
        print('model list:')
        print(model_list)

    # create a list of the indices

    index_mesh = np.array(np.meshgrid(*[range(n) for n in shape]))

    index_list =  index_mesh.T.reshape(n_models, n_axes)
    if verbose: 
        print('index list:')
        print(index_list)

    return n_axes, shape, n_models, mesh, model_list, index_list


