#!/bin/bash

# create a cloudy grid using the default assumptions

# synthesizer_dir="/Users/sw376/Dropbox/Research/data/synthesizer/"
synthesizer_dir="/research/astrodata/highz/synthesizer/" # apollo
machine="apollo"
incident_cloudy_model="blackbody"
c="/research/astro/flare/software/cloudy/"
cloudy_params="c17.03-limited" 

cd ..
python create_cloudy_input_grid_from_cloudy.py -synthesizer_data_dir $synthesizer_dir -machine $machine -incident_cloudy_model $incident_cloudy_model  -cloudy_params $cloudy_params  -cloudy_path $c

