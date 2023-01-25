#!/usr/bin/env python3

from config import config, bounds
import sys
import os

config = config() # initialise variables in python

# IF USING A DEVELOPMENT BRANCH OF COAST, ADD THE REPOSITORY TO PATH:
#sys.path.append('/home/h01/fred/NOTES/SET_UP_CONDA_ARTIFACTORY/COAsT')
#sys.path.append('/home/h01/fred/NOTES/SET_UP_CONDA_ARTIFACTORY/COAST_SCIPY')
#sys.path.append('/data/users/fred/SET_UP_CONDA_ARTIFACTORY/COAST_SCIPY')

print(config.coast_repo)
sys.path.append(config.coast_repo)

import coast
#import xarray as xr
#import numpy as np
#import datetime
#import pandas as pd

print('Modules loaded')
args = sys.argv

year = int(args[1])
month  = int(args[2])


fn_en4 = config.din_en4 + "EN.4.2.2.f.profiles.g10.%d%02d.nc"%(year, month)
print(f"Should add filename specification to config file")

# Make target dir if not present
print(os.popen("mkdir -p "+config.dout_en4).read())
fn_out = config.dout_en4 + config.region+"_processed_%d%02d.nc"%(year, month)
print (f"Load file: {fn_en4}")

profile = coast.Profile(config=config.fn_cfg_prof)
profile.read_en4(fn_en4, multiple=True)

# Restrict by region
b = bounds(config.region)
profile = profile.subset_indices_lonlat_box(lonbounds=b.lonbounds, latbounds=b.latbounds)

# Process data: apply QC flags
try:
    new_profile = profile.process_en4(remove_flagged_neighbours=True)
except:
    new_profile = profile.process_en4()
    print(f"You were supposed to update COAsT!!")
    break

new_profile.dataset.to_netcdf(fn_out)
print(f"Saved: {fn_out}")
