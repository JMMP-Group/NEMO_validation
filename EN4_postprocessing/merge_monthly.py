import xarray as xr
import sys
from dask.diagnostics import ProgressBar
from config import config

config = config() # initialise variables in python

args = sys.argv
model = args[1]  # MOD. Already loaded into config.dn_out directory path
month = int(args[2])  # month

# Merge over all available years: "????" are 4-digit year labels
ds_index = xr.open_mfdataset(config.dn_out +
                             'interpolated_profiles_p0_????%02d_????.nc'%(month),
                             combine='nested', concat_dim="id_dim", parallel=True)
ds_diff = xr.open_mfdataset(config.dn_out +
                            'profile_errors_p0_????%02d_????.nc'%(month),
                            combine='nested', concat_dim="id_dim", parallel=True)


with ProgressBar():
  ds_index.to_netcdf(config.dn_out+"%02d_PRO_INDEX.nc"%(month))
  ds_diff.to_netcdf(config.dn_out+"%02d_PRO_DIFF.nc"%(month))

print('File written to {config.dn_out+"%02d_PRO_INDEX.nc"%(month)}')
print('File written to {config.dn_out+"%02d_PRO_DIFF.nc"%(month)}')




