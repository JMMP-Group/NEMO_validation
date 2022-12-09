#!/usr/bin/env python3

from config import config, bounds
import sys
import os

config = config() # initialise variables in python

# IF USING A DEVELOPMENT BRANCH OF COAST, ADD THE REPOSITORY TO PATH:
#sys.path.append('/home/h01/fred/NOTES/SET_UP_CONDA_ARTIFACTORY/COAsT')
#sys.path.append('/home/h01/fred/NOTES/SET_UP_CONDA_ARTIFACTORY/COAST_SCIPY')
# sys.path.append('/data/users/fred/SET_UP_CONDA_ARTIFACTORY/COAST_SCIPY')

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

#fn_en4 = "/scratch/fred/EN4/EN.4.2.2.f.profiles.g10.*.nc"
#fn_en4 = "/scratch/fred/EN4/EN.4.2.2.f.profiles.g10.%d%02d.nc"%(year,month)
fn_en4 = config.din_en4 + "EN.4.2.2.f.profiles.g10.%d%02d.nc"%(year,month)

# Make target dir if not present
print(os.popen("mkdir -p "+config.dout_en4).read())
#fn_out = "/scratch/fred/EN4/SCIPY_processed_%d%02d.nc"%(year,month)
fn_out = config.dout_en4 + config.region+"_processed_%d%02d.nc"%(year,month)
print (fn_en4)
print (fn_out)
#fn_cfg_prof = "/data/users/fred/coast_demo/config/example_en4_profiles.json"
#fn_cfg_prof = "/home/users/jelt/GitHub/COAsT/config/example_en4_profiles.json"

profile = coast.Profile(config=config.fn_cfg_prof)
profile.read_en4(fn_en4, multiple=True)

# Restrict by region
b = bounds(config.region)
profile = profile.subset_indices_lonlat_box(lonbounds=b.lonbounds, latbounds=b.latbounds)

new_profile = profile.process_en4()

new_profile.dataset.to_netcdf(fn_out)
print(f"Saved: {fn_out}")
