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

from config import config
import sys

config = config() # initialise variables in python

# IF USING A DEVELOPMENT BRANCH OF COAST, ADD THE REPOSITORY TO PATH:
sys.path.append(config.coast_repo)

import coast
import xarray as xr
import numpy as np
#import datetime
import pandas as pd
import os
from coast import crps_util as cu
import numpy as np

import time

def reduce_resolution(profile_mod, profile_obs):
    """
    Take the nearest x,y,t grid point information from Profile object that came from the model.
    Identify profiles from common grid points.
    Average observation profiles over common profiles (model does not have sufficient spatial/temporal resolution to
    resolve their differences).
    Output reduced set of obs profiles corresponding to unique model space-time points
    Output corresponding model profiles.

    NB If the model profiles were first interpolated on a single obs profile and then onto a common reference profile, then an
    assumption is being made that the vertical distribution of the first obs profile (in a set to be averaged) is
    representative of the whole set.

    Assumes the z_dim dimension and depth coordinate is common to all depth varying variables.
    I.e. all profiles are on the same vertical grid.

    INPUTS:
    profile_mod (z_dim, id_dim): output from Profile.obs_operator() method with nearest_x, nearest_y, nearest_t
                                    variables. A COAsT.Profile object
    profile_obs (z_dim, id_dim): A COAsT.Profile object

    OUTPUTS:
    profile_mod_reduced (z_dim, id_dim_reduced): A COAsT.Profile object with unique (nearest_x, nearest_y, nearest_t)
    profile_obs_reduced (z_dim, id_dim_reduced): A COAsT.Profile object where data has been averaged.

    COAsT.Profile object
    """

    ds_mod = profile_mod.dataset  # provides nearest grid point indices and new coords
    ds = profile_obs.dataset  # provides the obs dataset, to be reduced
    var_modifier = ""  # optional suffix for output variable names

    # Find the unique triples of nearest model x,y,t indices (within the Profile obj). And the indices where they are found (in the parent Profile object)
    unique_mod_indices, unique_mod_indices_id = np.unique(np.array([[int(ds_mod.nearest_index_x[i].values),
        int(ds_mod.nearest_index_y[i].values),
        int(ds_mod.nearest_index_t[i].values)] for i in range(ds_mod.dims['id_dim'])]), axis=0, return_index=True) # find unique rows


    # Get a list of variables in this dataset
    vars_in = [items for items in ds.keys()]
    vars_out = [vv + "{0}".format(var_modifier) for vv in vars_in]

    # Get output dimensions and create 2D id, depth arrays
    n_id = len(unique_mod_indices)
    n_z = ds.dims['z_dim']

    # Create output dataset and fill in new coordinates
    depth_arr = ds_mod.depth.isel(id_dim=unique_mod_indices_id).values
    time_arr = ds_mod.time.isel(id_dim=unique_mod_indices_id).values
    lat_arr = ds_mod.latitude.isel(id_dim=unique_mod_indices_id).values
    lon_arr = ds_mod.longitude.isel(id_dim=unique_mod_indices_id).values

    ds_out = xr.Dataset(coords=dict(depth=(["id_dim", "z_dim"], depth_arr),
                                    time=(["id_dim"], time_arr),
                                    longitude=(["id_dim"], lon_arr),
                                    latitude=(["id_dim"], lat_arr)))

    # Loop over variables and create empty placeholders
    for vv in vars_out:
        if len(ds[vv].shape) == 2:
            ds_out[vv] = (["id_dim", "z_dim"], np.zeros((n_id, n_z)) * np.nan)
        elif len(ds[vv].shape) == 1:
            ds_out[vv] = (["id_dim"], np.zeros((n_id)) * np.nan)
        else:
            print(f"Panic! Did not expect vv:{vv}, shape:{ds[vv].shape}")

    # Grid_N is the count in each box
    ds_out["grid_N{0}".format(var_modifier)] = (["id_dim", "z_dim"], np.zeros((n_id, n_z)) * np.nan)

    # Loop over every unique profile (at model resolution)
    print(f"Averaging obs onto the space-time model resolution - slow")
    for ii in range(n_id):
        condition_x = np.isclose(ds_mod.nearest_index_x.values, unique_mod_indices[ii][0])
        condition_y = np.isclose(ds_mod.nearest_index_y.values, unique_mod_indices[ii][1])
        condition_t = np.isclose(ds_mod.nearest_index_t.values, unique_mod_indices[ii][2])
        condition_xy = np.logical_and(condition_x, condition_y)
        prof_ind = np.logical_and(condition_t, condition_xy)  # all the profile indices for THIS nearest model point
        # Check the profile indices include unique_mod_indices_id
        if prof_ind[unique_mod_indices_id[ii]] == False:
            print(f"Panic. Arrays have become unordered. Expected ii:{ii}, unique_mod_indices_id[ii]=True")
        # NB Thoughts on speed up options:
        # prof_ind is an array of np.bool(np.size(input id_dim))
        # unique_mod_indices_id are integers
        # Such that np.argmin(prof_ind>0) = unique_mod_indices_id, but unique... does not capture repeats.
        # Using unique_mod_indices_id could speed up the algorithm but .isel(id_dim=prof_ind).mean(dim="id_dim") fails
        # as isel collapses the dimension for one value

        for vv in range(len(vars_in)):
            vv_in = vars_in[vv]
            vv_out = vars_out[vv]

            # Check datatype isn't string
            if ds[vv_in].dtype in ["S1", "S2"]:
                continue
            ds_out[vv_out][ii] = ds[vv_in].isel(id_dim=prof_ind).mean(dim="id_dim")
            # if len(ds[vv_in].shape) == 2:
            #    ds_out[vv_out][ii] = ds[vv_in].isel(id_dim=prof_ind).mean(dim="id_dim")
            # elif len(ds[vv_in].shape) == 1:
            #    ds_out[vv_out][ii] = ds[vv_in].isel(id_dim=prof_ind).mean(dim="id_dim")

            # Store N in own variable
            ds_out["grid_N{0}".format(var_modifier)][ii] = np.sum(prof_ind)

    print(f"Averaging obs onto the space-time model resolution - DONE")

    # Subselect model profiles - assumes model can be interpolated onto the depth bins of the
    #  first obs depth profiles, as representative of the depth bins for the averaged profiles
    ds_mod_decimated = ds_mod.isel(id_dim=unique_mod_indices_id)
    print(f"Subselect model profiles onto the space-time model resolution - DONE")

    return coast.Profile(dataset=ds_mod_decimated), coast.Profile(dataset=ds_out)

###########################################################

def surface_crps_process(gridded_mod_surf, prof_obs_surf):

    """
    gridded_mod_surf  xr.Dataset with temperature and salinity xr.Dataarrays
    prof_obs_surf  xr.Dataset with temperature and salinity xr.Dataarrays and latotide, longitude, time coords (or variables)

    """
    radius_list = [0, 8, 14, 20]  # evaluate CRPS over radii (km)

    var_list = ["temperature", "salinity"]
    n_id = prof_obs_surf.dims['id_dim']
    n_rad = len(radius_list)
    n_var = len(var_list)

    for v_count, var_str in enumerate(var_list):
        crps_vals = np.zeros((n_var, n_rad, n_id))*np.nan
        crps_points = np.zeros((n_var, n_rad, n_id), dtype=int)
        crps_land_flags = np.full((n_var, n_rad, n_id), True)
        for r_count, nh_radius in enumerate(radius_list):
            crps_vals[v_count, r_count,:], \
            crps_points[v_count, r_count,:], \
            crps_land_flags[v_count, r_count,:] = cu.crps_sonf_moving(gridded_mod_surf[var_str],
                            prof_obs_surf.longitude, prof_obs_surf.latitude, prof_obs_surf[var_str], prof_obs_surf.obs_time,
                            nh_radius,
                            'nearest')


    print(f"CRPS values: {crps_vals}")
    #print(f"Number of points used: {b}")
    #print(f"Land present?: {~land_flags}")
    #print(f"Average {var_str} CRPS (rad:{nh_radius}) where no land: {np.nanmean(crps_vals[~land_flags])}")

    # Add the crps metrics along new dimension
    for v_count, var_str in enumerate(var_list):
        prof_obs_surf[var_str+"_crps"] = (['radius','id_dim'], np.array(crps_vals[v_count,:,:]))
        prof_obs_surf[var_str+"_crps_pts"] = (['radius','id_dim'], np.array(crps_points[v_count,:,:]))
        prof_obs_surf[var_str+"_crps_land_flags"] = (['radius','id_dim'], np.array(crps_land_flags[v_count,:,:]))
    # Add coords to new dimension
    prof_obs_surf = prof_obs_surf.assign_coords({"radius": np.array(radius_list)})

    #print(prof_obs_surf.temperature_crps)
    print(prof_obs_surf)
    return prof_obs_surf


starttime =time.perf_counter()


args = sys.argv

exper = args[1]
startyear=int(args[2])
month=int(args[3])
endyear=int(args[4])
grid=args[5]
try:
	debug_flag = str(args[6])=="debug"
except: debug_flag = False

print('Modules loaded')

# Start and end dates for the analysis. The script will cut down model
# and EN4 data to be witin this range.
start_date = np.datetime64(str(startyear)+"-01-01")
end_date = np.datetime64(str(endyear)+"-01-01")
#end_date = np.datetime64(str(startyear)+"-02-01")


# Reference depths (in metres)
ref_depth = np.concatenate((np.arange(1,100,2), np.arange(100,300,5), np.arange(300, 1000, 50), np.arange(1000,4000,100)))

# Name of the run -- used mainly for naming output files
run_name='p0_%d%02d_%d'%(startyear,month,endyear)

# File paths (All)
fn_dom = "%s%s"%(config.dn_dom, grid)


# Say a month at a time
fn_dat = "%s%s%02d*T.nc"%(config.dn_dat, startyear, month)  # NB config.dn_dat contains $MOD/exper  ## NEEDS TO MOVE TO CONFIG
#fn_dat = "%scoast_example_nemo_subset_data.nc"%(config.dn_dat)  # NB config.dn_dat contains $MOD/exper
print(fn_dat)

dn_out = f"{config.dn_out}"

# Make them in case they are not there.
print(os.popen(f"mkdir -p {dn_out}").read())

fn_prof = config.dout_en4 + config.region+"_processed_%d%02d.nc"%(startyear,month)  # generated by pre_process_en4_monthly.py
fn_cfg_nemo = config.fn_cfg_nemo
fn_cfg_prof = config.fn_cfg_prof

# CREATE NEMO OBJECT and read in NEMO data. Extract latitude and longitude array
print('Reading model data..')
print(f"nemo = coast.Gridded({fn_dat}, {fn_dom}, multiple=True, config={fn_cfg_nemo})")
nemo = coast.Gridded(fn_dat, fn_dom, multiple=True, config=fn_cfg_nemo)
print(nemo.dataset)
lon = nemo.dataset.longitude.values.squeeze()
lat = nemo.dataset.latitude.values.squeeze()
print('NEMO object created')

# Extract time indices between start and end dates
t_ind = nemo.dataset.time.values>=start_date
nemo.dataset = nemo.dataset.isel(t_dim=t_ind)
t_ind = nemo.dataset.time.values<end_date

BEFORE = starttime
NOW = time.perf_counter()
ALLTIME = NOW-starttime
DT = NOW-BEFORE
print("THIS FAR A %s %s ",ALLTIME,DT)
from sys import exit


nemo.dataset = nemo.dataset.isel(t_dim=t_ind)

# Create a landmask array -- important for obs_operator. We can 
# get a landmask from bottom_level.
nemo.dataset["landmask"] = nemo.dataset.bottom_level == 0
nemo.dataset = nemo.dataset.rename({"depth_0": "depth"})
print('Landmask calculated')

# CREATE EN4 PROFILE OBJECT containing processed data. We just need to
# create a Profile object and place the data straight into its dataset
profile = coast.Profile(config=fn_cfg_prof)
try:
  profile.dataset = xr.open_dataset(fn_prof, chunks={'id_dim':10000})
  print('Profile object created')
  # Extract time indices between start and end dates for Profile data.
  t_ind = pd.to_datetime(profile.dataset.time.values) >= start_date
  profile.dataset = profile.dataset.isel(id_dim=t_ind)
  t_ind = pd.to_datetime(profile.dataset.time.values) < end_date
  profile.dataset = profile.dataset.isel(id_dim=t_ind)
  print('Profile object time subsetted')
  print(profile.dataset)
except:
  profile.dataset = xr.open_dataset(fn_prof, chunks={'profile':10000})
  profile.dataset = profile.dataset.rename({"profile":"id_dim"})
  print(f"Profiles generated with legacy COAsT version. New dims: (id_dim, z_dim)")
  print('Profile object created')
  # Extract time indices between start and end dates for Profile data.
  t_ind = pd.to_datetime(profile.dataset.time.values)>=start_date
  profile.dataset = profile.dataset.isel(id_dim=t_ind)
  t_ind = pd.to_datetime(profile.dataset.time.values)<end_date
  profile.dataset = profile.dataset.isel(id_dim=t_ind)
  print('Profile object time subsetted')
  print(profile.dataset)

# Extract only the variables that we want
nemo.dataset = nemo.dataset[["temperature","salinity","bathymetry","bottom_level","landmask"]]
profile.dataset = profile.dataset[['potential_temperature','practical_salinity','depth']]
profile.dataset = profile.dataset.rename({"potential_temperature":"temperature", "practical_salinity":"salinity"})

#nemo.dataset = nemo.dataset[["temperature","bathymetry","bottom_level","landmask"]]
#profile.dataset = profile.dataset[['potential_temperature','depth']]
#profile.dataset = profile.dataset.rename({"potential_temperature":"temperature"})

# Cut out a geographical box - to speed up obs_operator processing
profile = profile.subset_indices_lonlat_box(lonbounds = [-26, 17],
					    latbounds = [44, 65])

# Cut out a time window - to speed up obs_operator processing
profile = profile.time_slice( date0=start_date, date1=end_date )

## Extract profiles from model
##############################
BEFORE = NOW
NOW = time.perf_counter()
ALLTIME = NOW-starttime
DT = NOW-BEFORE
print("THIS FAR B %s %s ",ALLTIME,DT)
#exit()

### TEMPORARILY RESTRICT NUMBER OF PROFILES WHILE DEBUGGING
print("profile.dataset[:]",profile.dataset)
if debug_flag == True:
	print(f"debug_flag = T # restrict number of profiles")
	#profile.dataset = profile.dataset.isel(id_dim=range(0,10))
	profile.dataset = profile.dataset.isel(id_dim=range(0,3000))

# Interpolate model to obs using obs_operator(). This is slow.
print("Extracting profiles from model - is slow ")

model_profiles = profile.obs_operator(nemo)
BEFORE = NOW
NOW = time.perf_counter()
ALLTIME = NOW-starttime
DT = NOW-BEFORE
print("THIS FAR C.0 %s %s ",ALLTIME,DT)

# Throw away profiles where the interpolation distance is larger than 5km.
keep_indices = model_profiles.dataset.interp_dist <= 5  ## SHOULD MOVE PARAMTER TO CONFIG
model_profiles = model_profiles.isel(id_dim=keep_indices)
profile = profile.isel(id_dim=keep_indices)
if np.sum(~keep_indices.values)>0:
	print(f"Dropped {np.sum(~keep_indices.values)} profiles: too far in space")

# Throw away profile where the interpolation time is larger than 12h
keep_indices = np.abs(model_profiles.dataset.interp_lag) <= np.timedelta64(12, 'h')  ## SHOULD MOVE PARAMTER TO CONFIG
model_profiles = model_profiles.isel(id_dim=keep_indices)
profile = profile.isel(id_dim=keep_indices)
if np.sum(~keep_indices.values)>0:
	print(f"Dropped {np.sum(~keep_indices.values)} profiles: too far in time")


BEFORE = NOW
NOW = time.perf_counter()
ALLTIME = NOW-starttime
DT = NOW-BEFORE
print("THIS FAR C.1 %s %s ",ALLTIME,DT)

# Load the profiles
#profile.dataset.load()
print('Model interpolated to obs locations')

## Perform analysis
###################
analysis = coast.ProfileAnalysis()

BEFORE = NOW
NOW = time.perf_counter()
ALLTIME = NOW-starttime
DT = NOW-BEFORE
print("THIS FAR C.2 %s %s ",ALLTIME,DT)

# Interpolate model profiles onto observation depths
model_profiles_interp = analysis.interpolate_vertical(model_profiles, profile, interp_method="linear")
print('Model interpolated to obs depths')

BEFORE = NOW
NOW = time.perf_counter()
ALLTIME = NOW-starttime
DT = NOW-BEFORE
print("THIS FAR C.3 %s %s ",ALLTIME,DT)

# Vertical interpolation of model profiles to reference depths
model_profiles_interp_ref_full = analysis.interpolate_vertical(model_profiles_interp, ref_depth)
print('Model interpolated to ref depths')

BEFORE = NOW
NOW = time.perf_counter()
ALLTIME = NOW-starttime
DT = NOW-BEFORE
print("THIS FAR C.4 %s %s ",ALLTIME,DT)

# Interpolation of obs profiles to reference depths
profile_interp_ref_full = analysis.interpolate_vertical(profile, ref_depth)
print('Obs interpolated to reference depths')

BEFORE = NOW
NOW = time.perf_counter()
ALLTIME = NOW-starttime
DT = NOW-BEFORE
print("THIS FAR C.5 %s %s ",ALLTIME,DT)

# Average depth-interpolated obs profiles onto horizontal-space/time model grid
#model_profiles_interp_ref, profile_interp_ref = analysis.reduce_resolution(model_profiles_interp_ref_full, profile_interp_ref_full)
model_profiles_interp_ref, profile_interp_ref = reduce_resolution(model_profiles_interp_ref_full, profile_interp_ref_full)

BEFORE = NOW
NOW = time.perf_counter()
ALLTIME = NOW-starttime
DT = NOW-BEFORE
print("THIS FAR C.6 %s %s ",ALLTIME,DT)

# Difference between Model and Obs
differences = analysis.difference(profile_interp_ref, model_profiles_interp_ref)

#differences.dataset.load()
print('Calculated errors')
BEFORE = NOW
NOW = time.perf_counter()
ALLTIME = NOW-starttime
DT = NOW-BEFORE
print("THIS FAR C.7 %s %s ",ALLTIME,DT)

# Surface & Bottom averaging
surface_def = 5
model_profiles_surface = analysis.depth_means(model_profiles_interp_ref, [0, surface_def])

BEFORE = NOW
NOW = time.perf_counter()
ALLTIME = NOW-starttime
DT = NOW-BEFORE
print("THIS FAR C.8 %s %s ",ALLTIME,DT)
obs_profiles_surface   = analysis.depth_means(profile_interp_ref, [0, surface_def])
surface_errors         = analysis.difference(obs_profiles_surface, model_profiles_surface)

BEFORE = NOW
NOW = time.perf_counter()
ALLTIME = NOW-starttime
DT = NOW-BEFORE
print("THIS FAR C %s %s ",ALLTIME,DT)

print(surface_errors.dataset)
print(model_profiles_surface.dataset)
print(obs_profiles_surface.dataset)
BEFORE = NOW
NOW = time.perf_counter()
ALLTIME = NOW-starttime
DT = NOW-BEFORE
print("THIS FAR D %s %s ",ALLTIME,DT)
surface_data = xr.merge((surface_errors.dataset, model_profiles_surface.dataset, obs_profiles_surface.dataset),
			   compat='override')
# Try Mid water, aiming for 1500m centered say 1200,1700
model_profiles_mid = analysis.depth_means(model_profiles_interp_ref, [1200, 1700])
obs_profiles_mid   = analysis.depth_means(profile_interp_ref, [1200, 1700])
mid_errors         = analysis.difference(obs_profiles_mid, model_profiles_mid)

print(surface_errors.dataset)
print(model_profiles_surface.dataset)
print(obs_profiles_surface.dataset)
BEFORE = NOW
NOW = time.perf_counter()
ALLTIME = NOW-starttime
DT = NOW-BEFORE
print("THIS FAR E %s %s ",ALLTIME,DT)
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

print(surface_errors.dataset)
print(model_profiles_surface.dataset)
print(obs_profiles_surface.dataset)
BEFORE = NOW
NOW = time.perf_counter()
ALLTIME = NOW-starttime
DT = NOW-BEFORE
print("THIS FAR F %s %s ",ALLTIME,DT)
print('Bottom and surface data estimated')

# Write datasets to file
model_profiles.dataset.to_netcdf(dn_out+"extracted_profiles_{0}.nc".format(run_name))
model_profiles_interp_ref.dataset.to_netcdf(dn_out + "interpolated_profiles_{0}.nc".format(run_name))
profile_interp_ref.dataset.to_netcdf(dn_out + "interpolated_obs_{0}.nc".format(run_name))

differences.dataset.to_netcdf(dn_out+"profile_errors_{0}.nc".format(run_name))
surface_data.to_netcdf(dn_out+"surface_data_{0}.nc".format(run_name))
mid_data.to_netcdf(dn_out+"mid_data_{0}.nc".format(run_name))
bottom_data.to_netcdf(dn_out+"bottom_data_{0}.nc".format(run_name))

print('Analysis datasets written to file')

BEFORE = NOW
NOW = time.perf_counter()
ALLTIME = NOW-starttime
DT = NOW-BEFORE
print("THIS FAR G %s %s ",ALLTIME,DT)

if(0):


  print('CRPS analysis')

  # CRPS analysis of surface fields
  gridded_mod_surf = nemo.dataset.where(nemo.dataset.depth <= surface_def).mean(dim="z_dim")

  BEFORE = NOW
  NOW = time.perf_counter()
  ALLTIME = NOW-starttime
  DT = NOW-BEFORE
  print("THIS FAR H %s %s ",ALLTIME,DT)

  surface_data_crps = surface_crps_process(gridded_mod_surf, surface_data)
  surface_data_crps.to_netcdf(dn_out+"surface_crps_data_{0}.nc".format(run_name))

  BEFORE = NOW
  NOW = time.perf_counter()
  ALLTIME = NOW-starttime
  DT = NOW-BEFORE
  print("THIS FAR I %s %s ",ALLTIME,DT)
  print('CRPS Analysis done and datasets written to file')

