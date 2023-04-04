import xarray as xr
import matplotlib.pyplot as plt

import sys


from dask.diagnostics import ProgressBar



plt.figure()
for year in range(2004,2015):
 print (year)
 for month in range(1,12):
   print (month)
   ds_index = xr.open_mfdataset('/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/P0.0/analysisb/interpolated_profiles_p0_%04d%02d_????.nc'%(year,month),combine='nested',concat_dim="id_dim",parallel=True)
   plt.plot((year+month/12.,year+month/12.),(0,ds_index.id_dim.size),marker="o")
   #print(ds_index.id_dim.size)
#ds_diff = xr.open_mfdataset('/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/analysisb/profile_errors_p0_????%02d_????.nc'%(model,month),combine='nested',concat_dim="id_dim",parallel=True)

plt.savefig("HOWMANY_PROFILES.png")
