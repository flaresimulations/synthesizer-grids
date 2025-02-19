"""
Functions for creating submission scripts for specific machines.
"""


def create_slurm_job_script(
        number_of_models,
        partition,
        new_grid_name,
        cloudy_output_dir,
        cloudy_path,
        memory="4G"):

    """
    Create a generic slurm input script.
    """

    slurm_job_script = f"""#!/bin/bash
#SBATCH --job-name=run_cloudy      # Job name
#SBATCH --output=output/%A_%a.out  # Standard output log
#SBATCH --error=output/%A_%a.err   # Error log
#SBATCH --array=1-{number_of_models}  # Job array range
#SBATCH --ntasks=1                 # Number of tasks per job
#SBATCH --cpus-per-task=1          # CPU cores per task
#SBATCH --mem={memory}             # Memory per task
#SBATCH --partition={partition}    # Partition/queue name

# Run command
python run_cloudy.py \\
    --grid_name={new_grid_name} \\
    --cloudy_output_dir={cloudy_output_dir} \\
    --cloudy_path={cloudy_path} \\
    --index=${{SLURM_ARRAY_TASK_ID}}
"""
    return slurm_job_script


def artemis(
    photoionisation_n_models,
    partition,
    new_grid_name,
    cloudy_output_dir,
    cloudy_path,
    memory="4G"):

    """
    Submission script generator for artemis
    """

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

    # create job script
    slurm_job_script = create_slurm_job_script(
        photoionisation_n_models,
        partition,
        new_grid_name,
        cloudy_output_dir,
        cloudy_path,
        memory=memory)

    # save job script
    open(f"{new_grid_name}.slurm", "w").write(slurm_job_script)

