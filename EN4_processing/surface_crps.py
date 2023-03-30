#!/usr/bin/env python3
#
# This script will use CRPS analysis tools in COAsT for
# comparison between observed profile data to model data.
#
# I don't recommend applying this script to a whole dataset if the model
# run is long. I found it best to implement this script as part of a 
# parallelised routine using a slurm job array or similar. Each
# process will do a different time period (so start and end dates is
# all that needs to change between processes).
#
# To run: e.g.
# python surface_crps.py P0.0 2004 1 2005 CO7_EXACT_CFG_FILE.nc

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

###########################################################

def surface_crps_process(gridded_mod_surf, prof_obs_surf):

    """
    gridded_mod_surf  xr.Dataset with temperature and salinity xr.Dataarrays
    prof_obs_surf  xr.Dataset with temperature and salinity xr.Dataarrays and latitude, longitude, time coords (or variables)

    """
    radius_list = [0, 8, 14, 20]  # evaluate CRPS over radii (km)

    var_list = ["temperature", "salinity"]
    n_id = prof_obs_surf.dims['id_dim']
    n_rad = len(radius_list)
    n_var = len(var_list)

    crps_vals = np.zeros((n_var, n_rad, n_id))*np.nan
    crps_points = np.zeros((n_var, n_rad, n_id), dtype=int)
    crps_land_flags = np.full((n_var, n_rad, n_id), True)

    for v_count, var_str in enumerate(var_list):
        for r_count, nh_radius in enumerate(radius_list):
            print(f"{var_str}: **Radius**: {nh_radius} in {radius_list}")
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
nemo = coast.Gridded(fn_data=fn_dat, fn_domain=fn_dom, multiple=True, config=fn_cfg_nemo)
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


nemo.dataset = nemo.dataset.isel(t_dim=t_ind)

# Create a landmask array -- important for obs_operator. We can 
# get a landmask from bottom_level.
nemo.dataset["landmask"] = nemo.dataset.bottom_level == 0
nemo.dataset = nemo.dataset.rename({"depth_0": "depth"})
print('Landmask calculated')


# Extract only the variables that we want. NB NEMO temperature was mapped from potential temperature in the json file?
nemo.dataset = nemo.dataset[["temperature","salinity","bathymetry","bottom_level","landmask"]]

# Surface & Bottom averaging
surface_def = 5


if(1):
  #surface_data = xr.open_dataset(dn_out+"surface_data_{0}.nc".format(run_name), chunks={'id_dim': 10000})
  surface_data = xr.open_dataset(dn_out+"surface_data_{0}.nc".format(run_name))

  print('CRPS analysis')
  print(f"\n surface_data:\n {dn_out+'surface_data_{0}.nc'.format(run_name)}\n {surface_data}")

  # CRPS analysis of surface fields
  gridded_mod_surf = nemo.dataset.where(nemo.dataset.depth <= surface_def).mean(dim="z_dim")
  print(f"\n gridded_mod_surf:\n {gridded_mod_surf}")

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
  print(f'CRPS Analysis done and datasets written to file: {dn_out+"surface_crps_data_{0}.nc".format(run_name)}')
  print("Next merge and compute regional averages. E.g. merge_mean_surface_crps.py")
