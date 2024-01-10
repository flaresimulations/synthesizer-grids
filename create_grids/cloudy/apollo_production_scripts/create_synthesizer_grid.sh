#!/bin/bash
# synthesizer_dir="/Users/stephenwilkins/Dropbox/Research/data/synthesizer/" # SW's machine
synthesizer_dir="/research/astrodata/highz/synthesizer/" # apollo


cd ..
while IFS="" read -r p || [ -n "$p" ]
do
  arrIN=(${p// / })
  grid=${arrIN[0]}
  echo $grid
  python3 create_synthesizer_grid.py -grid_name $grid -synthesizer_data_dir $synthesizer_dir -line_calc_method linelist
done < apollo_production_scripts/grids.txt



