#!/bin/bash

# create a cloudy grid using the default assumptions

# synthesizer_dir="/Users/stephenwilkins/Dropbox/Research/data/synthesizer/"
synthesizer_dir="/research/astrodata/highz/synthesizer/" # apollo
machine="apollo"
c="/research/astro/flare/software/cloudy/"

cd ..
while IFS="" read -r p || [ -n "$p" ]
do
  arrIN=(${p// / })
  incident=${arrIN[0]}
  params=${arrIN[1]}
  printf '%s\n' "$sps"
  printf '%s\n' "$params"
  echo python3 create_cloudy_input_grid.py -synthesizer_data_dir $synthesizer_dir -machine $machine -incident_grid $incident  -cloudy_params $params  -cloudy_path $c
  python3 create_cloudy_input_grid.py -synthesizer_data_dir $synthesizer_dir -machine $machine -incident_grid $incident  -cloudy_params $params  -cloudy_path $c
done < apollo_production_scripts/incident.txt
