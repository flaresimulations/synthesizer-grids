
# Linelist generation

In our photoionisation modelling we need to provide cloudy with a list of lines which we're interested in. Here we use the outputs from several `cloudy` runs - including both stars and AGN and range of metallicities in both cases. Next we decide on a reference line and then gather lines that in at least one model are within some factor of the luminosity of that line. 

## Procedure

- extract stellar and AGN spectra for use by cloudy using `extract_spectra.ipynb`.
- run all the `cloudy` models. `cloudy` input files exist for both stellar and AGN sources, including both conventions for line labelling (naming), i.e. standard (air at >200nm, vacuum <200nm) or vacuum wavelengthd.
- run `create_linelist.ipynb` to generate the line list, changing the parameters as necessary.

