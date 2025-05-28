#!/usr/bin/env python3
from PythonEnvCfg.config import config
import sys
import numpy as np
import os

config = config() # initialise variables in python

# IF USING A DEVELOPMENT BRANCH OF COAST, ADD THE REPOSITORY TO PATH:
sys.path.append(config.coast_repo)

import xarray as xr
from dask.diagnostics import ProgressBar
import coast


# not used #

#def extract_season(ds, season=None):
#    # Extract season if needed
#    if season is not None:
#        season_array = coast.general_utils.determine_season(ds.time)
#        s_ind = season_array == season
#        ds = ds.isel(id_dim=s_ind)
#    return ds

##########

def merge_by_time(season="all"):
    """ merge all CRPS data into one file"""
    # set paths
    base_dir = config.dn_out + "profiles/"
    file_names = [f for f in os.listdir(base_dir) 
                  if f.startswith('surface_crps_data_')]
    
    ds_index = xr.open_mfdataset([config.dn_out + f for f in file_names],
                                 combine='nested', concat_dim="id_dim", 
                                 parallel=True)
    
    with ProgressBar():
      save_path = config.dn_out + "profiles/"
      ds_index.to_netcdf(save_path + "%03s_CRPS_merged.nc"%(season))
    print(f'File written to {config.dn_out+"%03s_CRPS_merged.nc"%(season)}')

def regional_means():
    """ mean over each region """

    # File paths (All)
    fn_dom_nemo = "%s%s"%(config.dn_dom, config.grid_nc)
    fn_cfg_nemo = config.fn_cfg_nemo
    fn_cfg_prof = config.fn_cfg_prof
    fn_analysis_crps = "%s%03s_CRPS_merged.nc"%(config.dn_out, season)
    fn_out = "%s%03s_mask_means_crps_daily.nc"%(config.dn_out, season)
    
    # get the CRPS data as a profile object
    crps = coast.Profile(config=fn_cfg_prof)
    crps.dataset = xr.open_dataset(fn_analysis_crps, chunks={'id_dim':10000})
    
    print('Doing regional analysis..')
    
    # load nemo lat/lon grid to define regions (as function of bathymetry)
    nemo = coast.Gridded(fn_domain=fn_dom_nemo, config=fn_cfg_nemo)
    
    bath = nemo.dataset.bathymetry.values.squeeze()
    lon = nemo.dataset.longitude.values.squeeze()
    lat = nemo.dataset.latitude.values.squeeze()
    
    # Define Regional Masks
    
    mm = coast.MaskMaker()
    masks_list = []
    # Add regional mask for whole domain
    masks_list.append(np.ones(lon.shape))  # 0
    masks_list.append(mm.region_def_nws_north_north_sea(lon, lat, bath))  # 1
    masks_list.append(mm.region_def_nws_outer_shelf(lon, lat, bath))  # 2
    masks_list.append(mm.region_def_nws_english_channel(lon, lat, bath))  # 3
    masks_list.append(mm.region_def_nws_norwegian_trench(lon, lat, bath))  # 4
    masks_list.append(mm.region_def_nws_kattegat(lon, lat, bath))  # 5
    masks_list.append(mm.region_def_nws_fsc(lon, lat, bath))  # 6
    #masks_list.append(mm.region_def_nws_shelf_break(lon,lat,bath))  # 7
    masks_list.append(mm.region_def_nws_south_north_sea(lon, lat, bath))  # 7
    masks_list.append(mm.region_def_nws_off_shelf(lon, lat, bath))  # 8
    masks_list.append(mm.region_def_nws_irish_sea(lon, lat, bath))  # 9
    
    masks_names = ['whole_domain', 'northern_north_sea','outer_shelf',
                   'eng_channel','nor_trench', 'kattegat', 'fsc', 
                   'southern_north_sea', 'off_shelf', 'irish_sea' ]
    print("Size of names is ", len(masks_names[:]))
    
    mask_xr = mm.make_mask_dataset(lon, lat, masks_list, masks_names)
    
    ## Perform analysis
    analysis = coast.ProfileAnalysis()
    mask_indices = analysis.determine_mask_indices(crps, mask_xr)
    
    # Add bathymetry data to profiles before mask averaging
    #differences.dataset['bathymetry'] = model_profiles_interp.dataset.bathymetry
    
    # Do mask averaging
    mask_means = analysis.mask_means(crps, mask_indices)
    # Rename averaged bathymetry variable. Drop duplicate
    #mask_means = mask_means.drop_vars("profile_mean_bathymetry").rename({"all_mean_bathymetry":"bathymetry"})
    
    # SAVE mask dataset to file
    #if season == "DJF": # only save once otherwise conflicts arise if writing same file simultantously
    #    mask_xr.to_netcdf(config.dn_out + 'mask_xr.nc')
    
    print('Regional means calculated.')
    
    print("mask means is", mask_means)
    
    # save mask dataset to file
    mask_means.to_netcdf(fn_out)
    print(f"Done. Saved mask_means to {fn_out}")

if __name__ == "__main__":
    merge_by_time()
    regional_means()
