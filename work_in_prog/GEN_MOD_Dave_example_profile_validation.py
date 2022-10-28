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

import sys
# IF USING A DEVELOPMENT BRANCH OF COAST, ADD THE REPOSITORY TO PATH:
#sys.path.append('/home/users/dbyrne/code/COAsT')
sys.path.append('/home/h01/fred/NOTES/SET_UP_CONDA_ARTIFACTORY/COAST_SCIPY')
import coast
import xarray as xr
import numpy as np
import datetime
import pandas as pd
import os

args = sys.argv

exper = args[1]
startyear=int(args[2])
endyear=int(args[3])
print('Modules loaded')

# Start and end dates for the analysis. The script will cut down model
# and EN4 data to be witin this range.
start_date = datetime.datetime(startyear,1,1)
end_date = datetime.datetime(endyear,1,1)

# Reference depths (in metres)
#ref_depth = np.concatenate((np.arange(1,100,2), np.arange(100,300,5), np.arange(300, 1000, 50)))
ref_depth = np.concatenate((np.arange(1,100,2), np.arange(100,300,5), np.arange(300, 1000, 50), np.arange(1000,4000,100)))
#ref_depth = np.arange(1,4000,10)

# Name of the run -- used mainly for naming output files
run_name='p0_%d_%d'%(startyear,endyear)

# File paths (All)
#fn_dom = "/gws/nopw/j04/jmmp_collab/CO9_AMM15/inputs/domains/CO7_EXACT_CFG_FILE.nc"
fn_dom = "/data/users/fred/CO7_EXACT_CFG_FILE.nc"
#fn_dat = "/gws/nopw/j04/jmmp/CO9_AMM15/outputs/{0}/daily/*.nc".format(run_name)
#fn_dat = "/gws/nopw/j04/jmmp_collab/AMM15/PORT/P0.0/*.nc"
#fn_dat = "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/DAILY/199[23]*T.nc"%(exper)
#fn_dat = "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/DAILY/%s*T.nc"%(exper,startyear)
fn_dat = "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/PPC3/%s*T.nc"%(exper,startyear)
print(fn_dat)
#dn_out = "/gws/nopw/j04/jmmp/CO9_AMM15/"
dn_out = "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/analysis/"%(exper)
# Make them in case they are not there.
print(os.popen("mkdir -p /scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/analysis"%(exper)).read())

#fn_prof = "/gws/nopw/j04/jmmp/CO9_AMM15/obs/processed.nc" # Processed eN4 file
fn_prof = "/scratch/fred/EN4/SCIPY_processed_1990-2020.nc"
#fn_cfg_nemo = "/home/users/dbyrne/enda/example_nemo_grid_t.json"
#fn_cfg_prof = "/home/users/dbyrne/enda/example_en4_profiles.json"

#fn_cfg_nemo = "/data/users/fred/coast_demo/config/example_nemo_grid_t_pot_pra.json"
fn_cfg_nemo = "/data/users/fred/coast_demo/config/example_nemo_grid_t.json"
fn_cfg_prof = "/data/users/fred/coast_demo/config/example_en4_profiles.json"
# CREATE NEMO OBJECT and read in NEMO data. Extract latitude and longitude array
print('Reading model data..')
nemo = coast.Gridded(fn_dat, fn_dom, multiple=True, config=fn_cfg_nemo)
print(nemo.dataset)
lon = nemo.dataset.longitude.values.squeeze()
lat = nemo.dataset.latitude.values.squeeze()
print('NEMO object created')

# Extract time indices between start and end dates
t_ind = pd.to_datetime(nemo.dataset.time.values)>=start_date
nemo.dataset = nemo.dataset.isel(t_dim=t_ind)
t_ind = pd.to_datetime(nemo.dataset.time.values)<end_date
nemo.dataset = nemo.dataset.isel(t_dim=t_ind)

# Create a landmask array -- important for obs_operator. We can 
# get a landmask from bottom_level.
nemo.dataset["landmask"] = nemo.dataset.bottom_level == 0
nemo.dataset = nemo.dataset.rename({"depth_0": "depth"})
print('Landmask calculated')

# CREATE EN4 PROFILE OBJECT containing processed data. We just need to
# create a Profile object and place the data straight into its dataset
profile = coast.Profile()
profile.dataset = xr.open_dataset(fn_prof, chunks={'profile':10000})
print('Profile object created')

# Extract time indices between start and end dates for Profile data.
t_ind = pd.to_datetime(profile.dataset.time.values)>=start_date
profile.dataset = profile.dataset.isel(profile=t_ind)
t_ind = pd.to_datetime(profile.dataset.time.values)<end_date
profile.dataset = profile.dataset.isel(profile=t_ind)

# Extract only the variables that we want
nemo.dataset = nemo.dataset[["temperature","salinity","bathymetry","bottom_level","landmask"]]
profile.dataset = profile.dataset[['potential_temperature','practical_salinity','depth']]
profile.dataset = profile.dataset.rename({"potential_temperature":"temperature", "practical_salinity":"salinity"})

# Interpolate model to obs using obs_operator()
model_profiles = profile.obs_operator(nemo)

# Throw away profiles where the interpolation distance is larger than 5km.
keep_indices = model_profiles.dataset.interp_dist <= 5
model_profiles = model_profiles.isel(profile=keep_indices)
profile = profile.isel(profile=keep_indices)

# Load the profiles
profile.dataset.load()
print('Model interpolated to obs locations')

# Vertical Interpolation of model profiles to obs depths
model_profiles_interp = model_profiles.interpolate_vertical(profile, interp_method="linear")
print('Model interpolated to obs depths')

# Vertical interpolation of model profiles to reference depths
model_profiles_interp = model_profiles_interp.interpolate_vertical(ref_depth)
print('Model interpolated to ref depths')

# Interpolation of obs profiles to reference depths
profile_interp = profile.interpolate_vertical(ref_depth)
print('Obs interpolated to reference depths')

# Difference between Model and Obs
differences = profile_interp.difference(model_profiles_interp)
differences.dataset.load()
print('Calculated errors')

# Surface & Bottom averaging
surface = 5
model_profiles_surface = model_profiles.depth_means([0, surface])
obs_profiles_surface = profile.depth_means([0, surface])
surface_errors = obs_profiles_surface.difference(model_profiles_surface)
print(surface_errors.dataset)
print(model_profiles_surface.dataset)
print(obs_profiles_surface.dataset)
surface_data = xr.merge((surface_errors.dataset, model_profiles_surface.dataset, obs_profiles_surface.dataset),
                           compat='override')
# Try Mid water, aiming for 1500m centered say 1200,1700
model_profiles_mid = model_profiles.depth_means([1200, 1700])
obs_profiles_mid = profile.depth_means([1200, 1700])
mid_errors = obs_profiles_mid.difference(model_profiles_mid)
print(mid_errors.dataset)
print(model_profiles_mid.dataset)
print(obs_profiles_mid.dataset)
mid_data = xr.merge((mid_errors.dataset, model_profiles_mid.dataset, obs_profiles_mid.dataset),
                           compat='override')

model_profiles_bottom = model_profiles.bottom_means([10, 30, 100],[100, 500, np.inf])
obs_bathymetry = model_profiles.dataset["bathymetry"].values
profile.dataset["bathymetry"] = (["profile"], obs_bathymetry)
obs_profiles_bottom = profile.bottom_means([10, 30, 100], [100, 500, np.inf])
bottom_errors = obs_profiles_bottom.difference(model_profiles_bottom)
bottom_data = xr.merge((bottom_errors.dataset, model_profiles_bottom.dataset, obs_profiles_bottom.dataset),
                          compat="override")
print('Bottom and surface data estimated')

# Write datasets to file
model_profiles.dataset.to_netcdf(dn_out+"extracted_profiles_{0}.nc".format(run_name))
model_profiles_interp.dataset.to_netcdf(dn_out + "interpolated_profiles_{0}.nc".format(run_name))
profile_interp.dataset.to_netcdf(dn_out + "interpolated_obs_{0}.nc".format(run_name))
differences.dataset.to_netcdf(dn_out+"profile_errors_{0}.nc".format(run_name))
surface_data.to_netcdf(dn_out+"surface_data_{0}.nc".format(run_name))
mid_data.to_netcdf(dn_out+"mid_data_{0}.nc".format(run_name))
bottom_data.to_netcdf(dn_out+"bottom_data_{0}.nc".format(run_name))
print('Analysis datasets writted to file')

print('Doing regional analysis..')
# Define Regional Masks
mm = coast.MaskMaker()
regional_masks = []
bath = nemo.dataset.bathymetry.values
regional_masks.append( np.ones(lon.shape) )
regional_masks.append( mm.region_def_nws_north_sea(lon,lat,bath))
regional_masks.append( mm.region_def_nws_outer_shelf(lon,lat,bath))
regional_masks.append( mm.region_def_nws_english_channel(lon,lat,bath))
regional_masks.append( mm.region_def_nws_norwegian_trench(lon,lat,bath))
regional_masks.append( mm.region_def_kattegat(lon,lat,bath))
regional_masks.append( mm.region_def_fsc(lon,lat,bath))
regional_masks.append( mm.region_def_south_north_sea(lon,lat,bath))
off_shelf = mm.region_def_off_shelf(lon, lat, bath)
off_shelf[regional_masks[3].astype(bool)] = 0
off_shelf[regional_masks[4].astype(bool)] = 0
regional_masks.append(off_shelf)
regional_masks.append( mm.region_def_irish_sea(lon,lat,bath))

region_names = ['whole_domain', 'north_sea','outer_shelf','eng_channel','nor_trench',
                'kattegat', 'fsc','southern_north_sea', 'irish_sea', 'off_shelf' ]

mask_list = mm.make_mask_dataset(lon, lat, regional_masks)
mask_indices = model_profiles_interp.determine_mask_indices(mask_list)

# Do mask averaging
mask_means = differences.mask_means(mask_indices)
print('Regional means calculated.')

# SAVE mask dataset to file
mask_means.to_netcdf(dn_out + 'mask_means_daily_{0}.nc'.format(run_name))
print('done')

