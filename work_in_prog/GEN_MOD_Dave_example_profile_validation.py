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
#
# To run: e.g.
# python GEN_MOD_Dave_example_profile_validation.py P0.0 2004 2005 
# or
# python GEN_MOD_Dave_example_profile_validation.py P0.0 2004 2005 debug
#
# last field "debug" restricts the number profiles to make manageable for debugging

import sys
# IF USING A DEVELOPMENT BRANCH OF COAST, ADD THE REPOSITORY TO PATH:
#sys.path.append('/home/users/dbyrne/code/COAsT')
#sys.path.append('/home/h01/fred/NOTES/SET_UP_CONDA_ARTIFACTORY/COAST_SCIPY')
sys.path.append('/home/users/jelt/GitHub/COAsT')
import coast
import xarray as xr
import numpy as np
#import datetime
import pandas as pd
import os



args = sys.argv

exper = args[1]
startyear=int(args[2])
endyear=int(args[3])
try:
	debug_flag = str(args[4])=="debug"
except: debug_flag = False
print('Modules loaded')

# Start and end dates for the analysis. The script will cut down model
# and EN4 data to be witin this range.
#start_date = datetime.datetime(startyear,1,1)
#end_date = datetime.datetime(endyear,1,1)
start_date = np.datetime64(str(startyear)+"-01-01")
end_date = np.datetime64(str(endyear)+"-01-01")

# Reference depths (in metres)
#ref_depth = np.concatenate((np.arange(1,100,2), np.arange(100,300,5), np.arange(300, 1000, 50)))
ref_depth = np.concatenate((np.arange(1,100,2), np.arange(100,300,5), np.arange(300, 1000, 50), np.arange(1000,4000,100)))
#ref_depth = np.arange(1,4000,10)

# Name of the run -- used mainly for naming output files
run_name='p0_%d_%d'%(startyear,endyear)

# File paths (All)
fn_dom = "/gws/nopw/j04/jmmp_collab/CO9_AMM15/inputs/domains/CO7_EXACT_CFG_FILE.nc"
#fn_dom = "/data/users/fred/CO7_EXACT_CFG_FILE.nc"
#fn_dat = "/gws/nopw/j04/jmmp/CO9_AMM15/outputs/{0}/daily/*.nc".format(run_name)
#fn_dat = "/gws/nopw/j04/jmmp_collab/AMM15/PORT/P0.0/*.nc"
#fn_dat = "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/DAILY/199[23]*T.nc"%(exper)
#fn_dat = "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/DAILY/%s*T.nc"%(exper,startyear)
fn_dat = "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/PPC3/%s*T.nc"%(exper,startyear)
fn_dat = "/home/users/deazer/SAMPLE_DATA_COAST_WORKFLOW/%s/DAILY/%s0*T.nc"%(exper,startyear)
print(fn_dat)
#dn_out = "/gws/nopw/j04/jmmp/CO9_AMM15/"
dn_out = "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/analysis/"%(exper)
dn_out = "/home/users/jelt/tmp/"
# Make them in case they are not there.
#print(os.popen("mkdir -p /scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/analysis"%(exper)).read())
print(os.popen("mkdir -p /home/users/jelt/tmp/").read())

#fn_prof = "/gws/nopw/j04/jmmp/CO9_AMM15/obs/processed.nc" # Processed eN4 file
fn_prof = "/scratch/fred/EN4/SCIPY_processed_1990-2020.nc"
fn_prof = "/home/users/deazer/SAMPLE_DATA_COAST_WORKFLOW/SCIPY_processed_1990-2020.nc"
#fn_cfg_nemo = "/home/users/dbyrne/enda/example_nemo_grid_t.json"
#fn_cfg_prof = "/home/users/dbyrne/enda/example_en4_profiles.json"

#fn_cfg_nemo = "/data/users/fred/coast_demo/config/example_nemo_grid_t_pot_pra.json"
fn_cfg_nemo = "/data/users/fred/coast_demo/config/example_nemo_grid_t.json"
fn_cfg_nemo = "/home/users/jelt/GitHub/COAsT/config/example_nemo_grid_t.json"
fn_cfg_prof = "/data/users/fred/coast_demo/config/example_en4_profiles.json"
fn_cfg_prof = "/home/users/jelt/GitHub/COAsT/config/example_en4_profiles.json"
# CREATE NEMO OBJECT and read in NEMO data. Extract latitude and longitude array
print('Reading model data..')
nemo = coast.Gridded(fn_dat, fn_dom, multiple=True, config=fn_cfg_nemo)
print(nemo.dataset)
lon = nemo.dataset.longitude.values.squeeze()
lat = nemo.dataset.latitude.values.squeeze()
print('NEMO object created')

# Extract time indices between start and end dates
t_ind = nemo.dataset.time.values>=start_date
nemo.dataset = nemo.dataset.isel(t_dim=t_ind)
t_ind = nemo.dataset.time.values<end_date
nemo.dataset = nemo.dataset.isel(t_dim=t_ind)

# Create a landmask array -- important for obs_operator. We can 
# get a landmask from bottom_level.
nemo.dataset["landmask"] = nemo.dataset.bottom_level == 0
nemo.dataset = nemo.dataset.rename({"depth_0": "depth"})
print('Landmask calculated')

# CREATE EN4 PROFILE OBJECT containing processed data. We just need to
# create a Profile object and place the data straight into its dataset
profile = coast.Profile(config=fn_cfg_prof)
profile.dataset = xr.open_dataset(fn_prof, chunks={'profile':10000})
profile.dataset = profile.dataset.rename({"profile":"id_dim"})
print('Profile object created')

print(profile.dataset)

# Extract time indices between start and end dates for Profile data.
t_ind = pd.to_datetime(profile.dataset.time.values)>=start_date
profile.dataset = profile.dataset.isel(id_dim=t_ind)
t_ind = pd.to_datetime(profile.dataset.time.values)<end_date
profile.dataset = profile.dataset.isel(id_dim=t_ind)

# Extract only the variables that we want
nemo.dataset = nemo.dataset[["temperature","salinity","bathymetry","bottom_level","landmask"]]
profile.dataset = profile.dataset[['potential_temperature','practical_salinity','depth']]
profile.dataset = profile.dataset.rename({"potential_temperature":"temperature", "practical_salinity":"salinity"})

# Cut out a geographical box - to speed up obs_operator processing
profile = profile.subset_indices_lonlat_box(lonbounds = [-26, 17],
					    latbounds = [44, 65])

# Cut out a time window - to speed up obs_operator processing
profile = profile.time_slice( date0=start_date, date1=end_date )

## Extract profiles from model
##############################

### TEMPORARILY RESTRICT NUMBER OF PROFILES WHILE DEBUGGING
if debug_flag == True:
	print(f"debug_flag = T # restrict number of profiles")
	profile.dataset = profile.dataset.isel(id_dim=range(0,10))

# Interpolate model to obs using obs_operator(). This is slow.
print("Extracting profiles from model - is slow ")
model_profiles = profile.obs_operator(nemo)

# Throw away profiles where the interpolation distance is larger than 5km.
keep_indices = model_profiles.dataset.interp_dist <= 5
model_profiles = model_profiles.isel(id_dim=keep_indices)
profile = profile.isel(id_dim=keep_indices)

# Load the profiles
#profile.dataset.load()
print('Model interpolated to obs locations')

## Perform analysis
###################
analysis = coast.ProfileAnalysis()


# Interpolate model profiles onto observation depths
model_profiles_interp = analysis.interpolate_vertical(model_profiles, profile, interp_method="linear")
print('Model interpolated to obs depths')

# Vertical interpolation of model profiles to reference depths
model_profiles_interp_ref = analysis.interpolate_vertical(model_profiles_interp, ref_depth)
print('Model interpolated to ref depths')

# Interpolation of obs profiles to reference depths
profile_interp_ref = analysis.interpolate_vertical(profile, ref_depth)
print('Obs interpolated to reference depths')

# Difference between Model and Obs
differences = analysis.difference(profile_interp_ref, model_profiles_interp_ref)

#differences.dataset.load()
print('Calculated errors')

# Surface & Bottom averaging
surface_def = 5
model_profiles_surface = analysis.depth_means(model_profiles_interp_ref, [0, surface_def])
obs_profiles_surface   = analysis.depth_means(profile_interp_ref, [0, surface_def])
surface_errors         = analysis.difference(obs_profiles_surface, model_profiles_surface)
print(surface_errors.dataset)
print(model_profiles_surface.dataset)
print(obs_profiles_surface.dataset)
surface_data = xr.merge((surface_errors.dataset, model_profiles_surface.dataset, obs_profiles_surface.dataset),
			   compat='override')
# Try Mid water, aiming for 1500m centered say 1200,1700
model_profiles_mid = analysis.depth_means(model_profiles_interp_ref, [1200, 1700])
obs_profiles_mid   = analysis.depth_means(profile_interp_ref, [1200, 1700])
mid_errors         = analysis.difference(obs_profiles_mid, model_profiles_mid)
print(mid_errors.dataset)
print(model_profiles_mid.dataset)
print(obs_profiles_mid.dataset)
mid_data = xr.merge((mid_errors.dataset, model_profiles_mid.dataset, obs_profiles_mid.dataset),
			   compat='override')


bottom_height = [10, 30, 100]  # Average over bottom heights of 10m, 30m and 100m for...
bottom_thresh = [100, 500, np.inf]  # ...bathymetry depths less than 100m, 100-500m and 500-infinite
model_profiles_bottom = analysis.bottom_means(model_profiles_interp_ref, bottom_height, bottom_thresh)
# create a bathymetry variable to take advantage of 
profile_interp_ref.dataset["bathymetry"] = (["id_dim"], model_profiles_interp_ref.dataset["bathymetry"].values)
obs_profiles_bottom = analysis.bottom_means(profile_interp_ref, bottom_height, bottom_thresh)
bottom_errors = analysis.difference( obs_profiles_bottom, model_profiles_bottom)


bottom_data = xr.merge((bottom_errors.dataset, model_profiles_bottom.dataset, obs_profiles_bottom.dataset),
			  compat="override")
print('Bottom and surface data estimated')

# Write datasets to file
model_profiles.dataset.to_netcdf(dn_out+"extracted_profiles_{0}.nc".format(run_name))
model_profiles_interp_ref.dataset.to_netcdf(dn_out + "interpolated_profiles_{0}.nc".format(run_name))
profile_interp_ref.dataset.to_netcdf(dn_out + "interpolated_obs_{0}.nc".format(run_name))
differences.dataset.to_netcdf(dn_out+"profile_errors_{0}.nc".format(run_name))
surface_data.to_netcdf(dn_out+"surface_data_{0}.nc".format(run_name))
mid_data.to_netcdf(dn_out+"mid_data_{0}.nc".format(run_name))
bottom_data.to_netcdf(dn_out+"bottom_data_{0}.nc".format(run_name))
print('Analysis datasets writted to file')



## Create a child class of MaskMaker in order to create/prototype a new region
class MaskMaker_new(coast.MaskMaker):
    @classmethod
    def region_def_fsc(cls, longitude, latitude, bath):
        """
        Regional definition for Faroe Shetland Channel (Northwest European Shelf)
        Longitude, latitude and bath should be 2D arrays corresponding to model
        coordinates and bathymetry. Bath should be positive with depth.
        """
        vertices_lon = [-7.13, -9.72, -6.37, -0.45, -4.53 ]
        vertices_lat = [62.17, 60.6, 59.07, 61.945, 62.51 ]
 
        mask = cls.fill_polygon_by_lonlat(np.zeros(longitude.shape), longitude, latitude, vertices_lon, vertices_lat)
        mask = mask * (bath > 200) * (bath > 0) * (~np.isnan(bath))
        return mask



# Define Regional Masks
#mm = coast.MaskMaker()
mm = MaskMaker_new()
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
off_shelf[regional_masks[3].astype(bool)] = 0  # excludes english channel (wasn't in anyway...)
off_shelf[regional_masks[4].astype(bool)] = 0  # exludes norwegian trench
regional_masks.append(off_shelf)
regional_masks.append( mm.region_def_irish_sea(lon,lat,bath))

region_names = ['whole_domain', 'north_sea','outer_shelf','eng_channel','nor_trench',
                'kattegat', 'fsc','southern_north_sea', 'off_shelf', 'irish_sea' ]

mask_list = mm.make_mask_dataset(lon, lat, regional_masks)
mask_indices = analysis.determine_mask_indices(model_profiles_interp_ref, mask_list)

# Mask plotting
if(0):
	import matplotlib.pyplot as plt
	for count in range(len(region_names)):
		plt.pcolormesh( mask_list.longitude, mask_list.latitude, mask_list.mask.isel(dim_mask=count)) 
		plt.contour(lon,lat,bath, [10,200], colors=["w","w"])
		plt.title(region_names[count])
		plt.savefig(f"mask_{region_names[count]}.png")


# Do mask averaging
mask_means = analysis.mask_means(differences, mask_indices)
print('Regional means calculated.')

# SAVE mask dataset to file
mask_means.to_netcdf(dn_out + 'mask_means_daily_{0}.nc'.format(run_name))
print('done')


