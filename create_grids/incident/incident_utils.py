"""
Generic incident grid creation functions module. 
"""

import os
import sys
import h5py
import numpy as np


# import functions from grid_utils module
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from grid_utils import get_grid_properties_from_hdf5  # need this to work
