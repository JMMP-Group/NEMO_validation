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
# run is long. I found it best to implement this script as part of a
# parallelised routine using a slurm job array or similar. Each
# process will do a different time period (so start and end dates is
# all that needs to change between processes).

from config import config
import sys

config = config() # initialise variables in python
# IF USING A DEVELOPMENT BRANCH OF COAST, ADD THE REPOSITORY TO PATH:
#sys.path.append('/home/users/dbyrne/code/COAsT')
#sys.path.append('/data/users/fred/SET_UP_CONDA_ARTIFACTORY/COAST_SCIPY')
sys.path.append(config.coast_repo)

import coast
import xarray as xr
import numpy as np
import datetime
import pandas as pd

args = sys.argv

model = args[1]
month = int(args[2])

# File paths (All)
#fn_dom = "/data/users/fred/CO7_EXACT_CFG_FILE.nc"
fn_dom_nemo = "%s%s"%(config.dn_dom, config.grid_nc)
#fn_dat = "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/P0.0/DAILY/20050101*T.nc*"
fn_dat_nemo = "%s%s%02d*T.nc"%(config.dn_dat, "2004", month)  # NB config.dn_dat contains $MOD/exper. yyyy:str is just to get grid data from a valid file
print(fn_dat_nemo)
#fn_cfg_nemo = "/data/users/fred/coast_demo/config/example_nemo_grid_t.json"
fn_cfg_nemo = config.fn_cfg_nemo
#fn_cfg_prof = "/data/users/fred/coast_demo/config/example_en4_profiles.json"
fn_cfg_prof = config.fn_cfg_prof
fn_analysis_diff = "%s%02d_PRO_DIFF.nc"%(config.dn_out, month)
fn_analysis_index = "%s%02d_PRO_INDEX.nc"%(config.dn_out, month)
#fn_out = "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/analysisb/%02d_mask_means_daily.nc"%(model,month)
fn_out = "%s%02d_mask_means_daily.nc"%(config.dn_out, month)

differences = coast.Profile(config=fn_cfg_prof)  ## IF OLD COAST FALLS OVER HERE, REMOVE "config=fn_cfg_prof" HERE AND +5 lines
#fn_in = "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/analysis/ALL_PRO_DIFF.nc"%(model)
#differences.dataset = xr.open_dataset("/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/analysisb/%02d_PRO_DIFF.nc"%(model,month), chunks={'id_dim':10000})
differences.dataset = xr.open_dataset(fn_analysis_diff, chunks={'id_dim':10000})
#differences.dataset = xr.open_dataset("EXZ_ALL_PRO_DIFF.nc", chunks={'profile':10000})
model_profiles_interp = coast.Profile(config=fn_cfg_prof)
#model_profiles_interp.dataset = xr.open_dataset("/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/analysisb/%02d_PRO_INDEX.nc"%(model,month), chunks={'id_dim':10000})
model_profiles_interp.dataset = xr.open_dataset(fn_analysis_index, chunks={'id_dim':10000})
#model_profiles_interp.dataset = xr.open_dataset("EXZ_ALL_PRO_INDEX.nc", chunks={'profile':10000})




print('Doing regional analysis..')

# load nemo lat/lon grid to define regions (as function of bathymetry)
#nemo = coast.Gridded(fn_dat_nemo, fn_domain=fn_dom_nemo, multiple=True, config=fn_cfg_nemo)
nemo = coast.Gridded(fn_domain=fn_dom_nemo, config=fn_cfg_nemo)

bath = nemo.dataset.bathymetry.values.squeeze()
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

if(1):    # New COAsT, without FSC as a defined region
    # Define Regional Masks

    ## Create a child class of MaskMaker in order to create/prototype a new region
    class MaskMaker_new(coast.MaskMaker):
        @classmethod
        def region_def_fsc(cls, longitude, latitude, bath):
            """
            Regional definition for Faroe Shetland Channel (Northwest European Shelf)
            Longitude, latitude and bath should be 2D arrays corresponding to model
            coordinates and bathymetry. Bath should be positive with depth.
            """
            vertices_lon = [-7.13, -9.72, -6.37, -0.45, -4.53]
            vertices_lat = [62.17, 60.6, 59.07, 61.945, 62.51]

            mask = cls.fill_polygon_by_lonlat(np.zeros(longitude.shape), longitude, latitude, vertices_lon, vertices_lat)
            mask = mask * (bath > 200) * (bath > 0) * (~np.isnan(bath))
            return mask

    mm = MaskMaker_new()
    masks_list = []
    # Add regional mask for whole domain
    masks_list.append(np.ones(lon.shape))  # 0
    masks_list.append(mm.region_def_nws_north_sea(lon, lat, bath))  # 1
    masks_list.append(mm.region_def_nws_outer_shelf(lon, lat, bath))  # 2
    masks_list.append(mm.region_def_nws_english_channel(lon, lat, bath))  # 3
    masks_list.append(mm.region_def_nws_norwegian_trench(lon, lat, bath))  # 4
    masks_list.append(mm.region_def_kattegat(lon, lat, bath))  # 5
    masks_list.append(mm.region_def_fsc(lon, lat, bath))  # 6
    #masks_list.append(mm.region_def_nws_shelf_break(lon,lat,bath))  # 7
    masks_list.append(mm.region_def_south_north_sea(lon, lat, bath))  # 7
    masks_list.append(mm.region_def_off_shelf(lon, lat, bath))  # 8
    masks_list.append(mm.region_def_irish_sea(lon, lat, bath))  # 9

    masks_names = ['whole_domain', 'north_sea','outer_shelf','eng_channel','nor_trench',
                    'kattegat','fsc', 'southern_north_sea', 'off_shelf', 'irish_sea' ]
    print( " Size of names is ",len(masks_names[:]))

    mask_xr = mm.make_mask_dataset(lon, lat, masks_list, masks_names)

    ## Perform analysis
    analysis = coast.ProfileAnalysis()
    mask_indices = analysis.determine_mask_indices(model_profiles_interp, mask_xr)

    # Do mask averaging
    mask_means = analysis.mask_means(differences, mask_indices)

    # SAVE mask dataset to file
    mask_xr.to_netcdf(config.dn_out + 'mask_xr.nc')

else:  # OLD COAsT with FSC region defined internally
    break  # really should not be using old coast as the region definitions have changed
    # Define Regional Masks
    mm = coast.MaskMaker()
    regional_masks = []

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
    #regional_masks.append( mm.region_def_nws_shelf_break(lon,lat,bath)) #7
    #print( " Size of masks is ",len(regional_masks[:]))
    regional_masks.append( mm.region_def_south_north_sea(lon,lat,bath)) #7
    print( " Size of masks is ",len(regional_masks[:]))
    off_shelf = mm.region_def_off_shelf(lon, lat, bath) #8
    off_shelf[regional_masks[3].astype(bool)] = 0 #8
    off_shelf[regional_masks[4].astype(bool)] = 0 #8
    regional_masks.append(off_shelf) #8
    print( " Size of masks is ",len(regional_masks[:]))
    regional_masks.append( mm.region_def_irish_sea(lon,lat,bath)) # 9
    print( " Size of masks is ",len(regional_masks[:]))


    print( " Final Size of masks is ",len(regional_masks[:]))

    region_names = ['whole_domain', 'north_sea','outer_shelf','eng_channel','nor_trench',
                    'kattegat','fsc', 'southern_north_sea', 'off_shelf', 'irish_sea' ]
    print( " Size of names is ",len(region_names[:]))

    mask_list = mm.make_mask_dataset(lon, lat, regional_masks)
    print(" mask_list is ",mask_list)
    mask_indices = model_profiles_interp.determine_mask_indices(mask_list)
    print(" mask_indices is ",mask_indices)

    # Do mask averaging
    differences.dataset = differences.dataset.rename({"id_dim": "profile"})
    mask_means = differences.mask_means(mask_indices)
    print(f"Doing difference.mask_means analysis with legacy COAsT version.")

print('Regional means calculated.')

print("mask means is", mask_means)

# SAVE mask dataset to file
#mask_means.to_netcdf('/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/analysisb/%02d_mask_means_daily.nc'%(model,month))
mask_means.to_netcdf(fn_out)
#mask_means.to_netcdf('EXZ_ALL_mask_means_daily.nc')
print('done')

