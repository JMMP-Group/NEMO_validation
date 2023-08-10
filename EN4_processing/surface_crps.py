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
# python surface_crps.py 2004 1

from config import config
import sys

config = config() # initialise variables in python

# IF USING A DEVELOPMENT BRANCH OF COAST, ADD THE REPOSITORY TO PATH:
sys.path.append(config.coast_repo)

import coast
import xarray as xr
import numpy as np
import os
from coast import crps_util as cu
import numpy as np
import time
starttime = time.perf_counter()

def surface_crps_process(gridded_mod_surf, prof_obs_surf):

    """
    gridded_mod_surf  xr.Dataset with temperature and salinity xr.Dataarrays
    prof_obs_surf  xr.Dataset with temperature and salinity xr.Dataarrays and
    latitude, longitude, time coords (or variables)
    """

    radius_list = [0, 8, 14, 20]  # evaluate CRPS over radii (km)

    var_list = ["temperature", "salinity"]
    n_id = prof_obs_surf.dims['id_dim']
    n_rad = len(radius_list)
    n_var = len(var_list)

    # Add coords to new dimension
    prof_obs_surf = prof_obs_surf.assign_coords(
                                             {"radius": np.array(radius_list)})

    # TODO: remove loops by constructing crps input/output as xarray types
    for var_str in var_list:
        vals_list, pts_list, lf_list = [], [], []
        for nh_radius in radius_list:

            print(f"{var_str}: **Radius**: {nh_radius} in {radius_list}")
            
            # calculate crps metric
            crps_vals, crps_points, crps_land_flags = cu.crps_sonf_moving(
                            gridded_mod_surf[var_str],
                            prof_obs_surf.longitude,
                            prof_obs_surf.latitude,
                            prof_obs_surf[var_str],
                            prof_obs_surf.obs_time,
                            nh_radius,
                            'nearest')

            # assign crps coords/dims
            dims=['radius','id_dim']
            coords = {'radius':(['radius'], [nh_radius])}
            kwargs = dict(coords=coords, dims=dims)

            # convert numpy to xarray
            xr_vals = xr.DataArray(crps_vals[None,:], **kwargs)
            xr_points = xr.DataArray(crps_points[None,:], **kwargs)
            xr_land_flags = xr.DataArray(crps_land_flags[None,:], **kwargs)

            # append list of xarray DataArrays
            vals_list.append(xr_vals)
            pts_list.append(xr_points)
            lf_list.append(xr_land_flags)

        # add crps metrics to Dataset
        prof_obs_surf[var_str + "_crps"] = xr.concat(vals_list, 'radius')
        prof_obs_surf[var_str+"_crps_pts"] = xr.concat(pts_list, 'radius')
        prof_obs_surf[var_str+"_crps_land_flags"] = xr.concat(lf_list, 'radius')

    print(f"CRPS values: {crps_vals}")

    return prof_obs_surf

args = sys.argv
start_year=int(args[1])
month=int(args[2])

# Name of the run -- used mainly for naming output files
date = '%d%02d'%(start_year,month)

# open ungridded surface data + LOAD
surface_data = xr.load_dataset(config.dn_out+"surface_data_{0}.nc".format(date))
# temp

# open gridded model data + LOAD
mod_path = "gridded_model_surface_data_%d%02d.nc"%(start_year, month)
print (config.dn_out + mod_path)
gridded_mod_surf = xr.load_dataset(config.dn_out + mod_path)

print('CRPS analysis')
print(f"""\n surface_data:
          \n {config.dn_out+'surface_data_{0}.nc'.format(date)}
          \n {surface_data}""")

surface_crps = surface_crps_process(gridded_mod_surf, surface_data)

# assign case number
surface_crps = surface_crps.assign_attrs(dict(case=config.case))

# save
surface_crps.to_netcdf(config.dn_out+"surface_crps_data_{0}.nc".format(date))


print(f"""CRPS Analysis done and datasets written to file:
         {config.dn_out+"surface_crps_data_{0}.nc".format(date)}""")
print("""Next merge and compute regional averages.
      E.g. merge_mean_surface_crps.py""")
