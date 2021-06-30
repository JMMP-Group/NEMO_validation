'''
Routines for calculating average differences between model runs. These routines
will average in time (either over some defined time period or into seasonal
groups) and in depth, according to some user defined reference depths.

The two main routines are:
    
    1. time_averaged_differences()
    2. depth_averaged_difference_means()
    
The first routine opens data for two model runs and will average specified 
variables in time. This output is saved to a new netCDF file, which should be
provided as input to depth_averaged_difference_means(). This will take the
time averaged means and further average them into user defined depth bins.

Please see the docstrings for the individual routines for more information.

'''

import pandas as pd
import xarray as xr
import sys
sys.path.append('/work/dbyrne/code/COAsT')
import numpy as np
import os.path
import coast
from dask.diagnostics import ProgressBar

def write_ds_to_file(ds, fn, **kwargs):
    ''' 
    Simple netcdf writing routine which checks if file already exists first 
    '''
    if os.path.exists(fn):
        os.remove(fn)
    with ProgressBar():
        ds.to_netcdf(fn, **kwargs)
        
def time_averaged_differences(fn_nemo_data1, fn_nemo_data2, fn_nemo_domain,
                              fn_out, fn_out_seasonal, variables = None, 
                              start_date = None, end_date = None,
                              chunks = {'time_counter':100}):
    '''
    Calculates and writes time-averaged differences between two provided 
    datasets. Differences will be calculated as dataset2 - dataset1.
    
    *NOTE: Ensure that both model configurations use the same grid (and
    NEMO domain file). Additionally, they should also lie over the same time
    period and frequency. No interpolation or regridding will be done by this
    routine. If model runs have the same time frequency and domain, but
    differing time periods, start_date and end_date can be specified over common
    time periods.
    
    INPUTS:
        fn_nemo_data1 (str)     : Full path to nemo data file 1
        fn_nemo_data2 (str)     : Full path to nemo data file 2
        fn_nemo_domain (str)    : Full path to nemo domain file
        fn_out (str)            : Full path to output file for all-period averaging
        fn_out_seasonal (str)   : Full path to output file for seasonal averaging
        variables (list of str) : List or array of COAsT variable names to include
                                  If unspecified, all time dependent variables
                                  will be included.
        start_date (datetime)   : Start date for analysis
        end_date (datetime)     : End date for analysis
        
    OUTPUTS:
            Two output files will be written by this routine. The first will contain
            differences averaged over the entire model runs (fn_out). The second wil
            contain climatological (seasonal) means.
    '''

    print('Loading NEMO data', flush=True)
    nemo1 = coast.NEMO(fn_nemo_data1, fn_nemo_domain, multiple=True, chunks=chunks)
    nemo2 = coast.NEMO(fn_nemo_data2, fn_nemo_domain, multiple=True, chunks=chunks)
    
    nemo1 = nemo1.dataset
    nemo2 = nemo2.dataset
    
    nemo1 = nemo1[variables]
    nemo2 = nemo2[variables]
    
    print('Subsetting time', flush=True)
    if (start_date is not None) and (end_date is not None):
        nemo1_time = pd.to_datetime(nemo1.time.values)
        t_ind = np.logical_and(nemo1_time>=start_date, nemo1_time<=end_date)
        nemo1 = nemo1.isel(t_dim=t_ind)
        
        nemo2_time = pd.to_datetime(nemo2.time.values)
        t_ind = np.logical_and(nemo2_time>=start_date, nemo2_time<=end_date)
        nemo2 = nemo2.isel(t_dim=t_ind)
    
    print('Calculating Differences', flush=True)
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
                            latitude = (['y_dim', 'x_dim'], nemo1.latitude),
                            depth_0 = (['z_dim','y_dim', 'x_dim'], nemo1.depth_0)))
    
    ds_out_season = xr.Dataset(coords = dict(
                            longitude = (['y_dim', 'x_dim'], nemo1.longitude),
                            latitude = (['y_dim', 'x_dim'], nemo1.latitude),
                            season = (['season'], mean_diff_season.season),
                            depth_0 = (['z_dim','y_dim', 'x_dim'], nemo1.depth_0)))
    
    print('Populating Output Dataset', flush = True)
    for vv in variables:
        ds_out[vv+'_mean_diff'] = mean_diff[vv]
        ds_out[vv+'_mean_abs_diff'] = mean_abs_diff[vv]
        ds_out[vv+'_max_pos_diff'] = max_pos_diff[vv]
        ds_out[vv+'_min_neg_diff'] = min_neg_diff[vv]
        ds_out[vv+'_std_diff'] = std_diff[vv]
        
        ds_out_season[vv+'_mean_diff'] = mean_diff_season[vv]
        ds_out_season[vv+'_mean_abs_diff'] = mean_abs_diff_season[vv]
        ds_out_season[vv+'_max_pos_diff'] = max_pos_diff_season[vv]
        ds_out_season[vv+'_min_neg_diff'] = min_neg_diff_season[vv]
        ds_out_season[vv+'_std_diff'] = std_diff_season[vv]
        
    ds_out['start_date'] = start_date
    ds_out['end_date'] = end_date
    ds_out_season['start_date'] = start_date
    ds_out_season['end_date'] = end_date
        
    print('Writing to file', flush = True)
    write_ds_to_file(ds_out, fn_out)
    write_ds_to_file(ds_out_season, fn_out_seasonal)
    
def depth_averaged_difference_means(fn, fn_seasonal, depth_bins,
                                    fn_out, fn_out_seasonal):
    
    ds = xr.open_dataset(fn, chunks={'z_dim':10})
    ds_season = xr.open_dataset(fn_seasonal, chunks={'z_dim':51})
    
    ny = ds.dims['y_dim']
    nx = ds.dims['x_dim']
    nz = len(depth_bins)
    nbins = nz-1
    
    bin_widths = [depth_bins[ii+1] - depth_bins[ii] for ii in np.arange(0,nbins)]
    bin_mids = [depth_bins[ii] + .5*bin_widths[ii] for ii in np.arange(0,nbins)]
    
    D = ds.depth_0.values
    
    # Average none seasonal data
    variables = list(ds.keys())
    ds_out = xr.Dataset(coords = dict(
                            longitude = (['y_dim', 'x_dim'], ds.longitude),
                            latitude = (['y_dim', 'x_dim'], ds.latitude),
                            bin_mids = (['depth_bin'], bin_mids),
                            bin_widths = (['depth_bin'], bin_widths)))
    for vv in variables:
        ds_out[vv] = (['depth_bin', 'y_dim', 'x_dim'], np.zeros((nbins, ny, nx)) * np.nan)
        
    for dd in np.arange(1,nbins):
        print(dd)
        D_tmp = np.logical_and(D<depth_bins[dd], D>=depth_bins[dd+1])
        
        for vv in variables:
            v_tmp = ds[vv].values
            if len(v_tmp.shape)==3:
               v_tmp[D_tmp] = np.nan
               ds_out[vv][dd-1] = np.nanmean(v_tmp, axis=0)
            
    write_ds_to_file(ds_out, fn_out)
            
    # Average seasonal data
    variables = list(ds_season.keys())
    n_season = ds_season.dims['season']
    ds_out_season = xr.Dataset(coords = dict(
                            longitude = (['y_dim', 'x_dim'], ds_season.longitude),
                            latitude = (['y_dim', 'x_dim'], ds_season.latitude),
                            bin_mids = (['depth_bin'], bin_mids),
                            bin_widths = (['depth_bin'], bin_widths),
                            season = (['season'], ds_season.season)))
    for vv in variables:
        ds_out_season[vv] = (['season', 'depth_bin', 'y_dim', 'x_dim'], np.zeros((n_season, nbins, ny, nx)) * np.nan)
        
    for dd in np.arange(1,nbins):
        print(dd)
        D_tmp = np.logical_and(D<depth_bins[dd], D>=depth_bins[dd+1])
        
        for vv in variables:
            v_tmp = ds_season[vv].values
            if len(v_tmp.shape)==4:
                v_tmp[:,D_tmp] = np.nan
                ds_out_season[vv][:,dd-1] = np.nanmean(v_tmp, axis=1)
            
    write_ds_to_file(ds_out_season, fn_out_seasonal)
    
    return
