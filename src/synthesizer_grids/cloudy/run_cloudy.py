"""
"""

import os
import yaml
import argparse

if __name__ == "__main__":

    # Initial parser
    parser = argparse.ArgumentParser(
        description="Run the cloudy models",)

    # Add additional parameters which are specific to this script

    parser.add_argument(
        "--cloudy-output-dir",
        type=str,
        required=False,
        default=None)
    
    parser.add_argument(
        "--grid-name",
        type=str,
        required=False,
        default=None)

    # Path to cloudy directory (not the executable; this is assumed to
    # {cloudy}/{cloudy_version}/source/cloudy.exe)
    parser.add_argument(
        "--cloudy-executable-path",
        type=str,
        required=False,
        default=None)

    # Path to cloudy directory (not the executable; this is assumed to
    # {cloudy}/{cloudy_version}/source/cloudy.exe)
    parser.add_argument("--index", type=int, required=True)

    # Parse arguments
    args = parser.parse_args()

    # Shorthand for cloudy_output_dir
    output_directory = f"{args.cloudy_output_dir}/{args.grid_name}"
    incident_index = args.index

    # Open the file containing the details of the photoionisation model
    parameter_file = f"{output_directory}/grid_parameters.yaml"
    with open(parameter_file, "r") as file:
        parameters = yaml.safe_load(file)

    # set CLOUDY_DATA_PATH environment variable
    os.environ['CLOUDY_DATA_PATH'] = (
        f"{args.cloudy_executable_path}/{parameters['cloudy_version']}/data/:./")

    # change directory to the output directory
    os.chdir(f"{output_directory}/{incident_index}")

    # Loop over each photoionisation model
    for i in range(parameters['photoionisation_n_models']):

        input_file = f"{output_directory}/{incident_index}/{i}.in"
        cloudy_executable = (
            f"{args.cloudy_executable_path}/{parameters['cloudy_version']}"
            "/source/cloudy.exe")

        # run the cloudy job
        command = f'{cloudy_executable} -r {i}'
        print(command)
        os.system(command)
