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
month = int(args[2])


fn_dom = "/data/users/fred/CO7_EXACT_CFG_FILE.nc"
fn_dat = "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/P0.0/DAILY/20050101*T.nc*"
print(fn_dat)
fn_cfg_nemo = "/data/users/fred/coast_demo/config/example_nemo_grid_t.json"
fn_cfg_prof = "/data/users/fred/coast_demo/config/example_en4_profiles.json"

differences = coast.Profile()
fn_in = "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/analysis/ALL_PRO_DIFF.nc"%(model)
differences.dataset = xr.open_dataset("/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/analysisb/%02d_PRO_DIFF.nc"%(model,month), chunks={'id_dim':10000})
differences.dataset = differences.dataset.rename({"id_dim": "profile"})
#differences.dataset = xr.open_dataset("EXZ_ALL_PRO_DIFF.nc", chunks={'profile':10000})
model_profiles_interp = coast.Profile()
model_profiles_interp.dataset = xr.open_dataset("/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/analysisb/%02d_PRO_INDEX.nc"%(model,month), chunks={'id_dim':10000})
#model_profiles_interp.dataset = xr.open_dataset("EXZ_ALL_PRO_INDEX.nc", chunks={'profile':10000})




print('Doing regional analysis..')
# Define Regional Masks
mm = coast.MaskMaker()
regional_masks = []
nemo = coast.Gridded(fn_dat, fn_dom, multiple=True, config=fn_cfg_nemo)

bath = nemo.dataset.bathymetry.values
lon = nemo.dataset.longitude.values.squeeze()
lat = nemo.dataset.latitude.values.squeeze()

#!regional_masks.append( np.ones(lon.shape) )
#!regional_masks.append( mm.region_def_nws_north_sea(lon,lat,bath))
#!regional_masks.append( mm.region_def_nws_outer_shelf(lon,lat,bath))
#!regional_masks.append( mm.region_def_nws_english_channel(lon,lat,bath))
#!regional_masks.append( mm.region_def_nws_norwegian_trench(lon,lat,bath))
#!regional_masks.append( mm.region_def_kattegat(lon,lat,bath))
#!regional_masks.append( mm.region_def_nws_fsc(lon,lat,bath))
#!regional_masks.append( mm.region_def_south_north_sea(lon,lat,bath))
#!
#!off_shelf = mm.region_def_off_shelf(lon, lat, bath)
#!off_shelf[regional_masks[3].astype(bool)] = 0
#!off_shelf[regional_masks[4].astype(bool)] = 0
#!regional_masks.append(off_shelf)
#!regional_masks.append( mm.region_def_irish_sea(lon,lat,bath))
#!
#!region_names = ['whole_domain', 'north_sea','outer_shelf','eng_channel','nor_trench',
#!                'kattegat', 'fsc', 'southern_north_sea', 'irish_sea', 'off_shelf' ]
#!!
regional_masks.append( np.ones(lon.shape) ) # 0
print( " Size of masks is ",len(regional_masks[:]))
regional_masks.append( mm.region_def_nws_north_sea(lon,lat,bath)) # 1
print( " Size of masks is ",len(regional_masks[:]))
regional_masks.append( mm.region_def_nws_outer_shelf(lon,lat,bath)) # 2
print( " Size of masks is ",len(regional_masks[:]))
regional_masks.append( mm.region_def_nws_english_channel(lon,lat,bath)) #3
print( " Size of masks is ",len(regional_masks[:]))
regional_masks.append( mm.region_def_nws_norwegian_trench(lon,lat,bath)) #4
print( " Size of masks is ",len(regional_masks[:]))
regional_masks.append( mm.region_def_kattegat(lon,lat,bath)) #5
print( " Size of masks is ",len(regional_masks[:]))
regional_masks.append( mm.region_def_nws_fsc(lon,lat,bath)) #6
print( " Size of masks is ",len(regional_masks[:]))
regional_masks.append( mm.region_def_nws_shelf_break(lon,lat,bath)) #7
print( " Size of masks is ",len(regional_masks[:]))
regional_masks.append( mm.region_def_south_north_sea(lon,lat,bath)) #8
print( " Size of masks is ",len(regional_masks[:]))
off_shelf = mm.region_def_off_shelf(lon, lat, bath) #9
off_shelf[regional_masks[3].astype(bool)] = 0 #9
off_shelf[regional_masks[4].astype(bool)] = 0 #9
regional_masks.append(off_shelf) #10
print( " Size of masks is ",len(regional_masks[:]))
regional_masks.append( mm.region_def_irish_sea(lon,lat,bath)) # 10
print( " Size of masks is ",len(regional_masks[:]))


print( " Final Size of masks is ",len(regional_masks[:]))

region_names = ['whole_domain', 'north_sea','outer_shelf','eng_channel','nor_trench',
                'kattegat','fsc','shelf_break', 'southern_north_sea', 'off_shelf', 'irish_sea' ]
print( " Size of names is ",len(region_names[:]))

mask_list = mm.make_mask_dataset(lon, lat, regional_masks)
print(" mask_list is ",mask_list)
mask_indices = model_profiles_interp.determine_mask_indices(mask_list)
print(" mask_indices is ",mask_indices)

# Do mask averaging
mask_means = differences.mask_means(mask_indices)
print('Regional means calculated.')

print("mask means is",mask_means)

# SAVE mask dataset to file
mask_means.to_netcdf('/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/analysisb/%02d_mask_means_daily.nc'%(model,month))
#mask_means.to_netcdf('EXZ_ALL_mask_means_daily.nc')
print('done')

