#!/usr/bin/env python3
from PythonEnvCfg.config import config
import sys

config = config() # initialise variables in python

# IF USING A DEVELOPMENT BRANCH OF COAST, ADD THE REPOSITORY TO PATH:
#sys.path.append('/home/users/dbyrne/code/COAsT')
#sys.path.append('/home/h01/fred/NOTES/SET_UP_CONDA_ARTIFACTORY/COAST_SCIPY')
#sys.path.append('/data/users/fred/SET_UP_CONDA_ARTIFACTORY/COAST_NOV_2022_DEVELOP/COAsT/')
#sys.path.append(config.coast_repo)

import xarray as xr
from dask.diagnostics import ProgressBar

def extract_season(ds, season):
    """ reduce data to given season """
    ds = ds.where(ds.time.dt.season.compute() == season, drop=True)
    return ds

def _preprocess(ds_month):
    """ drop broadcasting of depth variable """
    # TODO: this should be done in GEN_MOD_Dave_example_profile_vali...
    ds_month["depth"] = ds_month.depth.isel(id_dim=0)
    return ds_month

args = sys.argv
season = str(args[1])  # season: 'DJF', 'MAM', 'JJA', SON'

# Merge over all available years: "????" are 4-digit year labels
# interpolated profiles
ds_index = xr.open_mfdataset(config.dn_out + 
                             "profiles/interpolated_profiles_*.nc",
                             combine='nested', concat_dim="id_dim",
                             parallel=True, preprocess=_preprocess)
print (ds_index)
ds_index = extract_season(ds_index, season)
print (ds_index)
print (kfjs)

# profile bias
ds_diff = xr.open_mfdataset(config.dn_out +
                            'profiles/profile_errors_*.nc',
                            combine='nested', concat_dim="id_dim",
                            parallel=True, preprocess=_preprocess)
ds_diff = extract_season(ds_diff, season)

# observational profiles
ds_obs = xr.open_mfdataset(config.dn_out +
                           'profiles/interpolated_obs_*.nc',
                           combine='nested', concat_dim="id_dim",
                           parallel=True, preprocess=_preprocess)
ds_obs = extract_season(ds_diff, season)

# save
with ProgressBar():
  ds_index.to_netcdf(config.dn_out+"profiles/"+"%03s_PRO_INDEX.nc"%(season))
  ds_diff.to_netcdf(config.dn_out+"profiles/"+"%03s_PRO_DIFF.nc"%(season))
  ds_obs.to_netcdf(config.dn_out+"profiles/"+"%03s_PRO_OBS.nc"%(season))

print(f'File written to {config.dn_out+"profiles/"+"%03s_PRO_INDEX.nc"%(season)}')
print(f'File written to {config.dn_out+"profiles/"+"%03s_PRO_DIFF.nc"%(season)}')
