#!/usr/bin/env python3
from PythonEnvCfg.config import config
import sys

config = config() # initialise variables in python

# IF USING A DEVELOPMENT BRANCH OF COAST, ADD THE REPOSITORY TO PATH:
#sys.path.append('/home/users/dbyrne/code/COAsT')
#sys.path.append('/home/h01/fred/NOTES/SET_UP_CONDA_ARTIFACTORY/COAST_SCIPY')
#sys.path.append('/data/users/fred/SET_UP_CONDA_ARTIFACTORY/COAST_NOV_2022_DEVELOP/COAsT/')
sys.path.append(config.coast_repo)

import xarray as xr
from dask.diagnostics import ProgressBar
from coast import general_utils


def extract_season(ds, season=None):
    # Extract season if needed
    if season is not None:
        season_array = general_utils.determine_season(ds.time)
        s_ind = season_array == season
        ds = ds.isel(id_dim=s_ind)
    return ds

args = sys.argv
model = args[1]  # MOD. Already loaded into config.dn_out directory path
season = str(args[2])  # season: 'DJF', 'MAM', 'JJA', SON'

# Merge over all available years: "????" are 4-digit year labels
ds_index = xr.open_mfdataset(config.dn_out +
                             'interpolated_profiles_p0_*.nc',
                             combine='nested', concat_dim="id_dim", parallel=True)
ds_index = extract_season(ds_index, season)

ds_diff = xr.open_mfdataset(config.dn_out +
                            'profile_errors_p0_*.nc',
                            combine='nested', concat_dim="id_dim", parallel=True)
ds_diff = extract_season(ds_diff, season)



with ProgressBar():
  ds_index.to_netcdf(config.dn_out+"%03s_PRO_INDEX.nc"%(season))
  ds_diff.to_netcdf(config.dn_out+"%03s_PRO_DIFF.nc"%(season))

print(f'File written to {config.dn_out+"%03s_PRO_INDEX.nc"%(season)}')
print(f'File written to {config.dn_out+"%03s_PRO_DIFF.nc"%(season)}')




