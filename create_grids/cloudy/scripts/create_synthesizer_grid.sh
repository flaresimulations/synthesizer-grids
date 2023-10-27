#!/bin/bash
# synthesizer_dir="/Users/stephenwilkins/Dropbox/Research/data/synthesizer/" # SW's machine
synthesizer_dir="/research/astrodata/highz/synthesizer/" # apollo

grid="bpass-2.2.1-bin_chabrier03-0.1,100.0_cloudy-c17.03"

cd ..
python3 create_synthesizer_grid.py -grid_name $grid -synthesizer_data_dir $synthesizer_dir -line_calc_method linelist
