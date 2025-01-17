# synthesizer-grids

This package includes scripts and tools for generating grids for use by the `synthesizer` synthetic observations package. There are two main modules: scripts for creating incident grids (grids containing the theoretical spectra produced by a model) and scripts for creating cloudy input files and then combining them to form a `synthesizer` grid.

## Installation

To use `synthesizer-grids` you'll first need to install it. To do this you can simply invoke the following at the root of the repo.

```sh
pip install .
```

This will install all submodules containing tools needed by the grid generation scripts.

## Incident

The first thing we need to do is create an incident grid. In some contexts this may be all you will need.

Unfortunately, since SPS codes adopt heterogenous approaches to everything the code for generating `synthesizer` incident grids is equally heterogenous. **However**, many helper tools are provided to standardise the process as much as possible. The result of running any of the generation scripts will be a standardised output file that can be further processed (e.g. by `cloudy`) or used within `synthesizer`.

The generation scripts can be found within `src/synthesizer-grids/incident` with sub-directories for each type of model (e.g. `sps` or `blackholes`). These generation scripts all take command line arguments. These commandline arguments can be seen by running the script with the help flag.

```sh
python <script>.py --help
```

Most models only take the standard arguments, but others can include other options. Here is an example `--help` output from the BC03 SPS model containing the default arguments.

```sh
usage: install_bc03.py [-h] --grid_dir GRID_DIR [--input_dir INPUT_DIR] [--download]

BC03 download and grid creation

options:
  -h, --help              show this help message and exit
  --grid_dir GRID_DIR     The grid file path
  --input_dir INPUT_DIR   The input file path
  --download, -d          Should the data be downloaded?
```

Note that all scripts will only download the data needed for the model if the download flag is set. This means the data need only be downloaded once for each model.

## Processing grids with CLOUDY

In the `cloudy/` directory are scripts to create `cloudy` input scripts and then process them to create `synthesizer` grids. There are two steps:

### Creating cloudy input grids

There are two approaches to creating cloudy input grids, using an **incident grid** (for example from an SPS model as above) or generating the grid using a cloudy shape command (e.g. blackbody). These scearnios are handled by two separate modules.

#### Using an incident grid

To create an input grid, we run `make_cloudy_input_grid.py`. This takes the following arguments:

`-synthesizer_data_dir`: The directory containing both the cloudy runs (`cloudy/`) and grids (`grids/`). **Note** by default new grids are placed in `synthesizer_data_dir/grids/dev/`.
`-m (--machine)`: the type of HPC on which you are running. If not provided, no submission script will be written.
`-incident_grid`: The name of the synthesizer incident grid on which you wish to run `cloudy`.
`-cloudy_params`: the YAML file containing `cloudy` parameters for modelling. See the examples for more details. These can be single values or arrays resulting in the creating of higher-dimensionality grids. Also included is the `cloudy` version.
`-cloudy_path`: the path to the `cloudy` directory. **Note** this is not the path to executable. This was done to allow multiple versions of cloudy to exist simultanously.

#### Using a cloudy shape command

As an alternative we can create a grid directly from using one of `cloudy`'s in-built shape commands (e.g. blackbody). To do this we need to provide a yaml file containing the name of the model (at the moment this is limited to `blackbody` and `agn`) with all the parameters, including any that are to be varied as a list.

#### The param YAML file

This contains the parameters that will be input into CLOUDY.

Users can provide arrays/lists of parameters they wish to iterate over. In this case, the code will create individual runs for each of the parameters specified, with the following naming convention:

{sps*grid}*{imf}_cloudy_{param}\_{value}

where {param} is the key in the param file of the parameter to be varied, and {value} is the value to be provided. If this is numeric, it will be converted to a string, and if it is negative, the '-' will be substituted for 'm'.

### Creating cloudy synthesizer grids

Next we need to combine the cloudy ouputs to create a new `synthesizer` grid. To do this we run `create_synthesizer_grid.py`. This takes the following arguments:

```bash
usage: create_synthesizer_grid.py [-h] --grid-dir GRID_DIR [--verbose] --cloudy-dir CLOUDY_DIR --incident-grid INCIDENT_GRID [--cloudy-grid CLOUDY_GRID] [--cloudy-params CLOUDY_PARAMS]
                                  [--cloudy-params-addition CLOUDY_PARAMS_ADDITION] [--include-spectra] [--replace] [--line-calc-method LINE_CALC_METHOD] [--machine MACHINE] [--norm-by-Q]

Create Synthesizer HDF5 grid files from cloudy outputs.

options:
  -h, --help            show this help message and exit
  --grid-dir GRID_DIR   The directory containing (or to contain) the grids.
  --verbose             Are we talking?
  --cloudy-dir CLOUDY_DIR
                        The directory containing each of the individual cloudy run outputs.
  --incident-grid INCIDENT_GRID
                        The name of the incident grid (the grid used to generate the cloudy runs).
  --cloudy-grid CLOUDY_GRID
                        The name of the new cloudy reprocessed grid (default to <incident_grid>_cloudy.hdf5).
  --cloudy-params CLOUDY_PARAMS
                        The path to the cloudy parameter file.
  --cloudy-params-addition CLOUDY_PARAMS_ADDITION
                        The path to the cloudy 'extra' parameter file.
  --include-spectra     Should the spectra be included in the grid?
  --replace             Should missing files be replaced?
  --line-calc-method LINE_CALC_METHOD
                        The method used to calculate the line fluxes (either 'lines' or 'linelist')
  --machine MACHINE     The machine used to run the cloudy runs (currently only supports apollo)
  --norm-by-Q, -Q       Should the grid be normalised by the specific ionising luminosity? (default is True)

```

### Scripts

The cloudy directory also includes a script folder `apollo_production_scripts/` for processing multiple grids automatically. Essentially there are two scripts, one for creating the cloudy inputs and a second for creating the ultimate synthesizer grids. These loop over the contents of `incident.txt` and `grids.txt` to get the input parameters.

These scripts will only work on Sussex's Apollo HPC system. Thus, if you want to use something similar please copy and adapt them as required. If you want them to form part of the package please make a new folder.

## Contributing

Contributions welcome. Please read the contributing guidelines [here](https://github.com/flaresimulations/synthesizer/blob/main/docs/CONTRIBUTING.md).
