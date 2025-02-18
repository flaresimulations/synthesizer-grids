"""
Create a grid of cloudy scripts based on combining an incident grid with a set
of cloudy parameters. Also creates a machine specific script.

Creates a folder for every incident grid point and then input files for every
individual photoionisation model.

python create_cloudy_input_grid.py --incident_grid_name=bpass-2.2.1-bin_chabrier03-0.1,300.0-ages:6,7,8 --grid-dir=/Users/sw376/Dropbox/Research/data/synthesizer/grids --cloudy_output_dir=/Users/sw376/Dropbox/Research/data/synthesizer/cloudy --cloudy_params=c23.01-sps

python create_cloudy_input_grid.py --incident_grid_name=bpass-2.2.1-bin_chabrier03-0.1,300.0-ages:6,7,8 --grid-dir=/its/home/sw376/astrodata/synthesizer/grids --cloudy_output_dir=//its/home/sw376/astrodata/synthesizer/cloudy --cloudy_params=c23.01-sps --cloudy_params_additional=test_suite/ionisation_parameter

"""

from pathlib import Path
import shutil
import numpy as np
import yaml
from synthesizer.abundances import (
    Abundances,
)
from synthesizer.grid import Grid
from synthesizer.photoionisation import cloudy17, cloudy23

# local modules
from utils import (
    apollo_submission_script,
    get_grid_props_cloudy,
    get_cloudy_params,
    # get_grid_properties,
)

from synthesizer_grids.parser import Parser
# from synthesizer_grids.grid_io import GridFile

if __name__ == "__main__":
    parser = Parser(description="Run a grid of incident cloudy models")

    # The name of the incident grid
    parser.add_argument("--incident_grid_name", type=str, required=True)

    # The path to the incident grid
    # parser.add_argument("--incident_grid_dir", type=str, required=True)
    
    # Path to directory where cloudy runs are stored and run
    parser.add_argument("--cloudy_output_dir", type=str, required=True)

    # The cloudy reference parameter set
    parser.add_argument(
        "--cloudy_params", 
        type=str, 
        required=True, 
        default="c23.01-sps"
    )

    # A second cloudy parameter set which supersedes the above. This is used
    # when considering variations on the default parameter set.
    parser.add_argument(
        "--cloudy_params_addition",
        type=str,
        required=False,
    )

    # Path to cloudy directory (not the executable; this is assumed to
    # {cloudy}/{cloudy_version}/source/cloudy.exe)
    parser.add_argument("--cloudy_path", type=str, required=False)

    # Machine (for submission script generation)
    parser.add_argument("--machine", type=str, required=False)

    
    # Parse arguments
    args = parser.parse_args()

    incident_grid_name = args.incident_grid_name
    incident_grid_dir = args.grid_dir
    cloudy_output_dir = args.cloudy_output_dir
    cloudy_params = args.cloudy_params
    cloudy_params_addition = args.cloudy_params_addition
    cloudy_path = args.cloudy_path
    machine = args.machine

    print(incident_grid_name)
    print(incident_grid_dir)
    print(cloudy_output_dir)
    print(cloudy_params)
    print(cloudy_params_addition)


    # Open the incident grid using synthesizer
    # this is necessary 
    incident_grid = Grid(
        incident_grid_name,
        grid_dir=incident_grid_dir,
        read_lines=False,
    )

    incident_axes = incident_grid.axes
    incident_axes_values = incident_grid._axes_values

    # Get properties of the incident grid
    (
        incident_n_axes,
        incident_shape,
        incident_n_models,
        incident_mesh,
        incident_model_list,
        incident_index_list,
    ) = get_grid_props_cloudy(incident_axes, incident_axes_values, verbose=True)

    # Load the cloudy parameters you are going to run
    fixed_photoionisation_params, variable_photoionisation_params = (
        get_cloudy_params(cloudy_params))

    # If an additional parameter set is provided supersede the default
    # parameters with these.
    if cloudy_params_addition:
        additional_photoionisation_fixed_params, additional_photoionisation_variable_params = get_cloudy_params(cloudy_params_addition)
        fixed_photoionisation_params = fixed_photoionisation_params | additional_photoionisation_fixed_params
        variable_photoionisation_params = variable_photoionisation_params | additional_photoionisation_variable_params

    print(fixed_photoionisation_params)
    print(variable_photoionisation_params)


    # doesn't matter about the ordering of these
    photoionisation_axes = list(variable_photoionisation_params.keys())
    
    # get properties of the photoionsation grid
    (
        photoionisation_n_axes,
        photoionisation_shape,
        photoionisation_n_models,
        photoionisation_mesh,
        photoionisation_model_list,
        photoionisation_index_list,
    ) = get_grid_props_cloudy(photoionisation_axes,
                            variable_photoionisation_params,
                            verbose=True)




    # Get name of new grid (concatenation of incident_grid and cloudy
    # parameter file)
    new_grid_name = f"{incident_grid_name}_cloudy-{cloudy_params}"

    # If an additional parameter set append this to the new grid name
    if args.cloudy_params_addition:
        # Ignore the directory part if it exists
        cloudy_params_addition_name = args.cloudy_params_addition.split("/")[
            -1
        ]
        # Append the new_grid_name with the additional name
        new_grid_name += "-" + cloudy_params_addition_name


    # define output directory
    output_directory = f'{cloudy_output_dir}/{new_grid_name}'

    # make output directories
    Path(output_directory).mkdir(parents=True, exist_ok=True)
    
    # Loop over all incident models, extract the spectral energy distribtions,
    # convert to the cloudy format, and save in the cloudy folder.
    lam = incident_grid.lam

    # Use cloudy version to initialise a synthesizer.cloudy object
    if fixed_photoionisation_params["cloudy_version"] == "c17.03":
        cloudy = cloudy17
    if fixed_photoionisation_params["cloudy_version"] == "c23.01":
        cloudy = cloudy23

    # Save all the parameters for use later
    parameters_to_save = fixed_photoionisation_params | variable_photoionisation_params | incident_axes_values
    parameters_to_save = parameters_to_save | {'incident_n_models': int(incident_n_models), 'photoionisation_n_models': int(photoionisation_n_models)}

    # Convert numpy arrays to lists for use in YAML file
    for key, value in parameters_to_save.items():
        if type(value) is np.ndarray:
            print(key)
            parameters_to_save[key] = value.tolist()

    # Save parameters_to_save dictionary to a YAML file
    with open(f"{output_directory}/grid_parameters.yaml", "w") as file:
        yaml.dump(parameters_to_save, file, default_flow_style=False)

    # If the ionisation_parameter_model is the reference model (i.e. not fixed)
    # save the value of the ionising photon luminosity at the reference grid point.
    if fixed_photoionisation_params["ionisation_parameter_model"] == "ref":
        
        # Initialize an empty list to store reference values
        reference_values = []

        # Iterate over the axes of the incident grid
        for k in incident_grid.axes:

            # We should throw an exception here if a reference value is not
            #  included. This will happen for grids with additional axes.

            # Append the corresponding reference value from fixed_params
            reference_values.append(fixed_photoionisation_params["reference_" + k])

        # Make dictionary to allow unpacking when getting the reference
        # grid point
        ref_dict = dict(zip(incident_grid.axes, reference_values))

        # Get the reference grid point using the adjusted reference values
        incident_ref_grid_point = incident_grid.get_grid_point(**ref_dict)

        # Get the reference ionising photon luminosity 
        reference_log10_specific_ionising_lum = incident_grid.log10_specific_ionising_lum["HI"][incident_ref_grid_point]

        print(reference_log10_specific_ionising_lum)


    for i, (incident_params_tuple, incident_index_tuple) in enumerate(zip(
        incident_model_list, incident_index_list)):

        # get dictionary of incident parameters
        incident_parameters = dict(zip(incident_axes, incident_params_tuple))

        # print(incident_parameters)

        # Create a directory for each incident grid point
        Path(f"output_directory/{i}").mkdir(parents=True, exist_ok=True)

        # Create a directory for system outputs 
        # NOTE: this may be machine specific
        Path(f"{output_directory}/{i}/output").mkdir(parents=True, exist_ok=True)

        # Extract incident SED from the grid.
        lnu = incident_grid.spectra['incident'][tuple(incident_index_tuple)]

        # Save incident SED
        np.save(f"{output_directory}/{i}/input.ssed", np.array([lam, lnu]))

        # Convert incided SED to cloudy format and save
        shape_commands = cloudy23.ShapeCommands.table_sed(
            'input',
            lam,
            lnu,
            output_dir=f"{output_directory}/{i}")

        # Copy linelist to this folder
        shutil.copyfile('linelist.dat', f'{output_directory}/{i}/linelist.dat')

        # If the reference ionionisation parameter model is used we want
        # to calculate the difference between the specific ionising photon
        # luminosity for this grid point and the reference grid point
        if fixed_photoionisation_params["ionisation_parameter_model"] == "ref":
            delta_log10_specific_ionising_luminosity = (incident_grid.log10_specific_ionising_lum["HI"][tuple(incident_index_tuple)] - reference_log10_specific_ionising_lum)


        # Create cloudy input file for every photoionisation grid point

        # Loop over all photoionisation models
        # TEST THIS FOR single models
        for j, (photoionisation_params_tuple, photoionisation_index_tuple) in enumerate(zip(
        photoionisation_model_list, photoionisation_index_list)):
            
            # Get a dictionary of the photoionisation parameters that are varying
            # for this grid point
            variable_photoionisation_parameters = dict(
                zip(photoionisation_axes, photoionisation_params_tuple))

            # Create dictionary of all photoionisation parameters
            photoionisation_parameters = (variable_photoionisation_parameters 
                                        | fixed_photoionisation_params)

            # Create dictionary of all parameters
            parameters = incident_parameters | photoionisation_parameters 

            # Create synthesizer.Abundance object
            abundances = Abundances(
                metallicity=float(parameters["metallicity"]),
                reference=parameters["reference_abundance"],
                alpha=parameters["alpha"],
                abundances=parameters["abundance_scalings"],
                depletion_model=parameters["depletion_model"],
                depletion_scale=parameters["depletion_scale"],
            )

            # Define the ionisation parameter and geometry

            # If the reference ionionisation parameter model is used
            if parameters["ionisation_parameter_model"] == "ref":

                # For spherical geometry the effective log10U is this
                if parameters["geometry"] == "spherical":
                    ionisation_parameter = 10**(
                        np.log10(parameters["reference_ionisation_parameter"])
                        + (1 / 3) * delta_log10_specific_ionising_luminosity
                    )

                # For spherical geometry the effective log10U is this.
                # The spherical-U model uses U as a direct input to cloudy
                # instead of calculating Q.
                elif parameters["geometry"] == "spherical-U":
                    ionisation_parameter = 10**(
                        np.log10(parameters["reference_ionisation_parameter"])
                        + (1 / 3) * delta_log10_specific_ionising_luminosity
                    )

                # For plane-parallel geometry the effective just scales with
                # log10_specific_ionising_luminosity
                elif parameters["geometry"] == "planeparallel":
                    ionisation_parameter = 10**(
                        np.log10(parameters["reference_ionisation_parameter"])
                        + delta_log10_specific_ionising_luminosity
                    )

                else:
                    # TODO: turn this into an exception
                    print(
                        f"""ERROR: do not understand geometry choice:
                        {parameters['geometry']}"""
                    )

            # If fixed ionisation parameter model is used:
            elif parameters["ionisation_parameter_model"] == "fixed":
                ionisation_parameter = parameters["ionisation_parameter"]

            else:
                # TODO: turn this into an exception
                print(
                    f"""ERROR: do not understand U model choice:
                    {parameters['ionisation_parameter_model']}"""
                )

            # Set ionisation_parameter to provide cloudy
            parameters["ionisation_parameter"] = ionisation_parameter

            # Include shape command to read SED
            shape_commands = [f'table SED "input.sed" \n']

            # Output dictionary of all parameters as a YAML file
            with open(f"{output_directory}/{i}/{j}.yaml", "w") as file:
                yaml.dump(parameters_to_save, file, default_flow_style=False)

            # Create cloudy input file
            cloudy.create_cloudy_input(
                str(j),
                shape_commands,
                abundances,
                output_dir=f"{output_directory}/{i}",
                **parameters,
            )

    if machine == 'artemis':

        # determine the partition to use:
        
        # short = 2 hours
        if photoionisation_n_models < 5:
            partition = 'short'
        
        # general = 8 hours
        elif photoionisation_n_models < 33:
            partition = 'general'
        
        # long = 8 days
        else:
            partition = 'long'

        slurm_job_script = f"""#!/bin/bash
#SBATCH --job-name=run_cloudy      # Job name
#SBATCH --output=output/%A_%a.out  # Standard output log (%A = job ID, %a = task ID)
#SBATCH --error=error/%A_%a.err    # Error log
#SBATCH --array=1-{int(photoionisation_n_models)}               # Job array range
#SBATCH --ntasks=1                 # Number of tasks per job
#SBATCH --cpus-per-task=1          # CPU cores per task
#SBATCH --mem=4G                   # Memory per task
#SBATCH --partition={partition}          # Partition/queue name

# Run command
python run_cloudy.py --grid_name={new_grid_name} --cloudy_output_dir={cloudy_output_dir} --cloudy_path={cloudy_path} --index=${{SLURM_ARRAY_TASK_ID}}
"""

        open(f"{new_grid_name}.slurm", "w").write(slurm_job_script)


