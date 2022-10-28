#!/usr/bin/env python3
#
# This script will use all of the available analysis tools in COAsT for
# comparison of observed profile data to model data.
#
# At the top, please check the start and end dates, reference depths to
# interpolate data onto the the run_name, which is essentially the name
# of the output files.
#
# I don't recommend applying this script to a whole dataset if the model
# run is long. I found it best to implement this script as part of a 
# parallelised routine using a slurm job array or similar. Each
# process will do a different time period (so start and end dates is
# all that needs to change between processes).

import sys
# IF USING A DEVELOPMENT BRANCH OF COAST, ADD THE REPOSITORY TO PATH:
#sys.path.append('/home/users/dbyrne/code/COAsT')
sys.path.append('/data/users/fred/SET_UP_CONDA_ARTIFACTORY/COAST_SCIPY')
import coast
import xarray as xr
import numpy as np
import datetime
import pandas as pd


args = sys.argv

model = args[1]


fn_dom = "/data/users/fred/CO7_EXACT_CFG_FILE.nc"
fn_dat = "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/u-cp277/DAILY/20050101*T.nc*"
print(fn_dat)
fn_cfg_nemo = "/data/users/fred/coast_demo/config/example_nemo_grid_t.json"
fn_cfg_prof = "/data/users/fred/coast_demo/config/example_en4_profiles.json"

differences = coast.Profile()
fn_in = "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/analysis/ALL_PRO_DIFF.nc"%(model)
differences.dataset = xr.open_dataset("/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/analysis/ALL_PRO_DIFF.nc"%(model), chunks={'profile':10000})
#differences.dataset = xr.open_dataset("EXZ_ALL_PRO_DIFF.nc", chunks={'profile':10000})
model_profiles_interp = coast.Profile()
model_profiles_interp.dataset = xr.open_dataset("/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/analysis/ALL_PRO_INDEX.nc"%(model), chunks={'profile':10000})
#model_profiles_interp.dataset = xr.open_dataset("EXZ_ALL_PRO_INDEX.nc", chunks={'profile':10000})




print('Doing regional analysis..')
# Define Regional Masks
mm = coast.MaskMaker()
regional_masks = []
nemo = coast.Gridded(fn_dat, fn_dom, multiple=True, config=fn_cfg_nemo)

bath = nemo.dataset.bathymetry.values
lon = nemo.dataset.longitude.values.squeeze()
lat = nemo.dataset.latitude.values.squeeze()

regional_masks.append( np.ones(lon.shape) )
regional_masks.append( mm.region_def_nws_north_sea(lon,lat,bath))
regional_masks.append( mm.region_def_nws_outer_shelf(lon,lat,bath))
regional_masks.append( mm.region_def_nws_english_channel(lon,lat,bath))
regional_masks.append( mm.region_def_nws_norwegian_trench(lon,lat,bath))
regional_masks.append( mm.region_def_kattegat(lon,lat,bath))
regional_masks.append( mm.region_def_nws_fsc(lon,lat,bath))
regional_masks.append( mm.region_def_south_north_sea(lon,lat,bath))
regional_masks.append( region_def_nws_fsc(lon,lat,bath))


off_shelf = mm.region_def_off_shelf(lon, lat, bath)
off_shelf[regional_masks[3].astype(bool)] = 0
off_shelf[regional_masks[4].astype(bool)] = 0
regional_masks.append(off_shelf)
regional_masks.append( mm.region_def_irish_sea(lon,lat,bath))

region_names = ['whole_domain', 'north_sea','outer_shelf','eng_channel','nor_trench',
                'kattegat', 'fsc', 'southern_north_sea', 'irish_sea', 'off_shelf' ]

mask_list = mm.make_mask_dataset(lon, lat, regional_masks)
mask_indices = model_profiles_interp.determine_mask_indices(mask_list)

# Do mask averaging
mask_means = differences.mask_means(mask_indices)
print('Regional means calculated.')

# SAVE mask dataset to file
mask_means.to_netcdf('/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/analysis/ALL_mask_means_daily.nc'%(model))
#mask_means.to_netcdf('EXZ_ALL_mask_means_daily.nc')
print('done')

