import xarray as xr
import xarray.ufuncs as uf
import sys
sys.path.append('/Users/dbyrne/code/COAsT/')
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import os.path
import coast.general_utils as gu
import coast

fn_nemo_data1 = "/Users/dbyrne/Projects/CO9_AMM15/data/nemo/20040103_shelftmb_grid_T.nc"
fn_nemo_data2 = "/Users/dbyrne/Projects/CO9_AMM15/data/nemo/20040102_shelftmb_grid_T.nc"
fn_nemo_domain = "/Users/dbyrne/Projects/CO9_AMM15/data/nemo/CO7_EXACT_CFG_FILE.nc"
fn_out = "/Users/dbyrne/Projects/CO9_AMM15/data/test.nc"
fn_out_seasonal = "/Users/dbyrne/Projects/CO9_AMM15/data/test_seasonal.nc"

variables = ['ssh']



def write_ds_to_file(ds, fn, **kwargs):
    ''' 
    Simple netcdf writing routine which checks if file already exists first 
    '''
    if os.path.exists(fn):
        os.remove(fn)
    ds.to_netcdf(fn, **kwargs)

nemo1 = coast.NEMO(fn_nemo_data1, fn_nemo_domain, chunks={'time_counter':100}).dataset
nemo2 = coast.NEMO(fn_nemo_data2, fn_nemo_domain, chunks={'time_counter':100}).dataset

nemo1 = nemo1[variables]
nemo2 = nemo2[variables]

diff = nemo2 - nemo1
diff['time'] = nemo1.time
diff = diff.set_coords(['time'])
abs_diff = abs(diff)
diff_neg = -diff

# Annual Differences
mean_diff = diff.mean(dim='t_dim', skipna = True)
std_diff = diff.std(dim='t_dim', skipna = True)
mean_abs_diff = abs_diff.mean(dim='t_dim', skipna = True)
max_pos_diff = diff.max(dim='t_dim', skipna = True)
min_neg_diff = diff_neg.max(dim='t_dim', skipna = True)

# Seasonal Differences
diff_season = diff.groupby('time.season')
abs_diff_season = abs_diff.groupby('time.season')
diff_neg_season = diff_neg.groupby('time.season')

mean_diff_season = diff_season.mean(dim='t_dim', skipna = True )
std_diff_season = diff_season.std(dim='t_dim', skipna = True)
mean_abs_diff_season = abs_diff_season.mean(dim='t_dim', skipna = True)
max_pos_diff_season = diff_season.max(dim='t_dim', skipna = True)
min_neg_diff_season = diff_neg_season.max(dim='t_dim', skipna = True)

ds_out = xr.Dataset(coords = dict(
                        longitude = (['y_dim', 'x_dim'], nemo1.longitude),
                        latitude = (['y_dim', 'x_dim'], nemo2.latitude)))

ds_out_season = xr.Dataset(coords = dict(
                        longitude = (['y_dim', 'x_dim'], nemo1.longitude),
                        latitude = (['y_dim', 'x_dim'], nemo2.latitude),
                        season = (['season'], mean_diff_season.season)))

for vv in variables:
    ds_out[vv+'_mean_diff'] = mean_diff[vv]
    ds_out[vv+'_mean_abs_diff'] = mean_abs_diff[vv]
    ds_out[vv+'max_pos_diff'] = max_pos_diff[vv]
    ds_out[vv+'min_neg_diff'] = min_neg_diff[vv]
    
    ds_out_season[vv+'_mean_diff'] = mean_diff_season[vv]
    ds_out_season[vv+'_mean_abs_diff'] = mean_abs_diff_season[vv]
    ds_out_season[vv+'max_pos_diff'] = max_pos_diff_season[vv]
    ds_out_season[vv+'min_neg_diff'] = min_neg_diff_season[vv]
    
write_ds_to_file(ds_out, fn_out)
write_ds_to_file(ds_out, fn_out_seasonal)
