#!/bin/bash

# create a cloudy grid using the default assumptions

# synthesizer_dir="/Users/sw376/Dropbox/Research/data/synthesizer/"
synthesizer_dir="/research/astrodata/highz/synthesizer/" # apollo
cloudy_dir=$synthesizer_dir/cloudy
grid_dir=$synthesizer_dir/grids/dev
machine="apollo"
incident_cloudy_model="agn"
c="/research/astro/flare/software/cloudy/"
cloudy_params="c17.03-blr" 

cd ..
python create_cloudy_input_grid_from_cloudy.py -grid_dir $grid_dir -cloudy_dir $cloudy_dir -machine $machine -incident_cloudy_model $incident_cloudy_model  -cloudy_params $cloudy_params  -cloudy_path $c

