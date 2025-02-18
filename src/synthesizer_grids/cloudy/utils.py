import numpy as np
import yaml


def get_cloudy_params(param_file="c23.01-sps", param_dir="params"):
    """
    Read cloudy parameters from a yaml parameter file

    Arguments:
        param_file (str)
            filename of the parameter file
        param_dir (str)
            directory containing the parameter file

    Returns:
        fixed_params (dict)
            dictionary of parameters that are fixed
        grid_params (dict)
            dictionary of parameters that vary on the grid
    """

    # open paramter file
    with open(f"{param_dir}/{param_file}.yaml", "r") as stream:
        try:
            params = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    grid_params = {}
    fixed_params = {}

    # Loop over parameters
    for k, v in params.items():
        # If parameter is a list store it in the grid_parameters dictionary
        # and convert to a numpy array
        if isinstance(v, list):
            grid_params[k] = list(map(float, v))

        # Otherwise store it in fixed_params dictionary
        else:
            fixed_params[k] = v

    return fixed_params, grid_params


def get_grid_props_cloudy(axes, axes_values, verbose=True):
    """
    Get the cloudy related properties of a grid.

    This returns some basic properties needed for processing an incident
    grid through CLOUDY.

    Args:
        axes (list)
            List of axes to vary.
        axes_values (dict)
            Dictionary of axis values.
        verbose (bool)
            Are we talking or not?

    Returns:
        n_axes (int)
            Number of axes.
        shape (list)
            Shape of the grid.
        n_models (int)
            Number of models.
        mesh (numpy.ndarray)
            Mesh of the grid.
        model_list (numpy.ndarray)
            List of models.
        index_list (numpy.ndarray)
            List of indices of each model within the grid.
    """
    # The grid axes
    if verbose:
        print(f"axes: {axes}")

    # Number of axes
    n_axes = len(axes)
    if verbose:
        print(f"number of axes: {n_axes}")

    # The shape of the grid (useful for creating outputs)
    shape = list([len(axes_values[axis]) for axis in axes])
    if verbose:
        print(f"shape: {shape}")

    # Determine number of models
    n_models = np.prod(shape)
    if verbose:
        print(f"number of models to run: {n_models}")

    # Create the mesh of the grid
    mesh = np.array(
        np.meshgrid(*[np.array(axes_values[axis]) for axis in axes])
    )

    # Create the list of the models
    model_list = mesh.T.reshape(n_models, n_axes)
    if verbose:
        print("model list:")
        print(model_list)

    # Create a list of the indices
    index_mesh = np.array(np.meshgrid(*[range(n) for n in shape]))
    index_list = index_mesh.T.reshape(n_models, n_axes)
    if verbose:
        print("index list:")
        print(index_list)

    return n_axes, shape, n_models, mesh, model_list, index_list




def apollo_submission_script(n, grid_data_dir, cloudy_path, cloudy_version):
    """
    Create an Apollo SGE submission script.

    Parameters
    ----------
    n : int
        Number of models to run, sets size of array job.
    synthesizer_data_dir : str
        where to write the submission script
    cloudy : str
        bash executable for CLOUDY

    Returns
    -------
    None
    """

    # cloudy executable
    cloudy = f"{cloudy_path}/{cloudy_version}/source/cloudy.exe"

    print(cloudy)

    # cloudy data dir
    cloudy_data_path = f'"{cloudy_path}/{cloudy_version}/data/:./"'

    apollo_job_script = f"""
######################################################################
# Options for the batch system
# These options are not executed by the script, but are instead read by the
# batch system before submitting the job. Each option is preceeded by '#$' to
# signify that it is for grid engine.
#
# All of these options are the same as flags you can pass to qsub on the
# command line and can be **overriden** on the command line. see man qsub for
# all the details
######################################################################
# -- The shell used to interpret this script
#$ -S /bin/bash
# -- Execute this job from the current working directory.
#$ -cwd
#$ -l h_vmem=4G
#$ -l m_mem_free=4G
# -- Job output to stderr will be merged into standard out. Remove this line if
# -- you want to have separate stderr and stdout log files
#$ -j y
#$ -o output/
# -- Send email when the job exits, is aborted or suspended
# #$ -m eas
# #$ -M YOUR_USERNAME@sussex.ac.uk
######################################################################
# Job Script
# Here we are writing in bash (as we set bash as our shell above). In here you
# should set up the environment for your program, copy around any data that
# needs to be copied, and then execute the program
######################################################################

# set cloudy data path
export CLOUDY_DATA_PATH={cloudy_data_path}

{cloudy} -r $SGE_TASK_ID
"""

    open(f"{grid_data_dir}/run_grid.job", "w").write(apollo_job_script)
    print(grid_data_dir)

    # define the qsub command
    qsub = f"qsub -t 1:{n} -q smp.q -pe openmp 2 run_grid.job "
    print(qsub)

    # create a script to execute the qsub command
    open(f"{grid_data_dir}/run.sh", "w").write(qsub)

    return
