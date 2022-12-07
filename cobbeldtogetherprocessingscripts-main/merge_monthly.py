import xarray as xr

import sys


from dask.diagnostics import ProgressBar

args = sys.argv
model = args[1]
month = int(args[2])

ds_index = xr.open_mfdataset('/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/analysisb/interpolated_profiles_p0_????%02d_????.nc'%(model,month),combine='nested',concat_dim="id_dim",parallel=True)
ds_diff = xr.open_mfdataset('/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/analysisb/profile_errors_p0_????%02d_????.nc'%(model,month),combine='nested',concat_dim="id_dim",parallel=True)


with ProgressBar():
  ds_index.to_netcdf("/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/analysisb/%02d_PRO_INDEX.nc"%(model,month))
  ds_diff.to_netcdf("/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/analysisb/%02d_PRO_DIFF.nc"%(model,month))


print('All Done')




