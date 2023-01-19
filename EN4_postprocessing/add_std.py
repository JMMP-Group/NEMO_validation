#!/usr/bin/env python3
#
# Add standard deviation to the profile diagnostics.
# E.g.:
# python add_std.py DJF

from config import config
import sys

config = config() # initialise variables in python
# IF USING A DEVELOPMENT BRANCH OF COAST, ADD THE REPOSITORY TO PATH:
sys.path.append(config.coast_repo)

import coast
import xarray as xr
import numpy as np

args = sys.argv

season = str(args[1])

file_in = "%s%03s_mask_means_daily.nc"%(config.dn_out, season)
file_out = "%s%03s_mask_means_daily2.nc"%(config.dn_out, season)
fn_cfg_prof = config.fn_cfg_prof

prof = coast.Profile(config=fn_cfg_prof)
prof.dataset = xr.open_dataset(file_in, chunks={'id_dim':10000})

print('Doing regional analysis..')
prof.dataset['profile_std_diff_temperature'] = np.sqrt(prof.dataset.profile_mean_square_diff_temperature - np.square(prof.dataset.profile_mean_diff_temperature))
prof.dataset['profile_std_diff_salinity'] = np.sqrt(prof.dataset.profile_mean_square_diff_salinity - np.square(prof.dataset.profile_mean_diff_salinity))


# SAVE mask dataset to file
prof.dataset.to_netcdf(fn_out)
print('done')

