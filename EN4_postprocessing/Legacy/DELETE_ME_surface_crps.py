#!/usr/bin/env python3
from config import config
import sys

config = config() # initialise variables in python

# IF USING A DEVELOPMENT BRANCH OF COAST, ADD THE REPOSITORY TO PATH:
sys.path.append(config.coast_repo)

import xarray as xr
from dask.diagnostics import ProgressBar
from coast import general_utils
import coast
from coast import crps_util as cu
import numpy as np

args = sys.argv
var_str = "temperature"
exper = args[1]  # not used?
startyear=int(args[2])
month=int(args[3])
endyear=int(args[4])
grid=args[5]

# Name of the run -- used mainly for naming output files
run_name='p0_%d%02d_%d'%(startyear,month,endyear)

# model data file paths
fn_dat = "%scoast_example_nemo_subset_data.nc"%(config.dn_dat)  # NB config.dn_dat contains $MOD/exper
fn_dom = "%s%s"%(config.dn_dom, grid)
fn_cfg_nemo = config.fn_cfg_nemo

# obs profile data file path (processed)
fn_surf = "{0}surface_data_{1}.nc".format(config.dn_out, run_name)  # surface_data_p0_200701_2008.nc



## Load NEMO surface data
nemo = coast.Gridded(fn_dat, fn_dom, multiple=True, config=fn_cfg_nemo).isel(z_dim=0)
mod_surf = nemo.dataset[var_str]  # DataArray with latitude, longitude, time coordinates

## Load EN4 profile (processed) surface data
## this file contains both model and profile data extracted for the surface layer. The raw variable (e.g. temperature)
## is derived from the observations not the model (because the model variable of the same name was overwritten in the process)
obs_surf = xr.open_dataset(fn_surf)

# Extract only the variables wanted
obs_surf = obs_surf[["latitude", "longitude", var_str, "obs_time"]]

"""
RECALL:

def crps_sonf_moving(
    mod_array: xr.DataArray,
    obs_lon: np.ndarray,
    obs_lat: np.ndarray,
    obs_var: np.ndarray,
    obs_time: np.ndarray,
    nh_radius: float,
    time_interp: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    
    Handles the calculation of single-observation neighbourhood forecast CRPS for a moving observation instrument.

    Differs from crps_sonf_fixed in that latitude and longitude are arrays of locations. Mod_array must contain
    dimensions x_dim, y_dim and t_dim and coordinates longitude, latitude, time.

    Args:
        mod_array (xr.DataArray):  DataArray from a Model Dataset.
        obs_lon (np.ndarray): Longitudes of fixed observation point.
        obs_lat (np.ndarray): Latitudes of fixed observation point.
        obs_var (np.ndarray): of floatArray of variable values, e.g time series.
        obs_time: (np.ndarray): of datetimeArray of times, corresponding to obs_var.
        nh_radius (float): Neighbourhood radius in km.
        time_interp (str): Type of time interpolation to use.

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray]: Array of CRPS values,
            Array containing the number of model points used for each CRPS value,
            Array of bools indicating where a model neighbourhood contained land.

"""
radius_list = [0, 8, 14, 20]
for nh_radius in radius_list:
    a,b,c = cu.crps_sonf_moving( mod_surf,
                         obs_surf.longitude, obs_surf.latitude, obs_surf[var_str], obs_surf.obs_time,
                         nh_radius,
                         'nearest')

    print(f"CRPS values: {a[~c]}")

    print(f"Number of points used: {b}")

    print(f"Land present?: {~c}")

    print(f"Average CRPS where no land: {np.nanmean(a[~c])}")




