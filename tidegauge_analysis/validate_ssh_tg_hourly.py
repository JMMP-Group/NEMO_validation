"""
module add jaspy
export CONDA_ENV="/home/users/jelt/miniconda3/envs/coast_dev"
source activate $CONDA_ENV
"""

# location of COAsT repo, if using a particular branch
coast_repo = "/home/users/jelt/GitHub/COAsT"

import sys
# IF USING A DEVELOPMENT BRANCH OF COAST, ADD THE REPOSITORY TO PATH:
sys.path.append(coast_repo)

import coast
from coast import general_utils as gu
from coast import plot_util as pu
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from datetime import datetime, timedelta
import pandas as pd
from dateutil.relativedelta import relativedelta
import matplotlib.dates as mdates
import utide as ut
import scipy.signal as signal
import os
import glob
from dask.diagnostics import ProgressBar

def write_ds_to_file(ds, fn, **kwargs):
    ''' 
    Simple netcdf writing routine which checks if file already exists first 
    '''
    if os.path.exists(fn):
        os.remove(fn)
    ds.to_netcdf(fn, **kwargs)

def compare_phase(g1, g2):
    g1 = np.array(g1)
    g2 = np.array(g2)
    r = (g1-g2)%360 - 360
    r[r<-180] = r[r<-180] + 360
    r[r>180] = r[r>180] - 360
    return r

def analyse_ssh_old(fn_ext, fn_out, thresholds = np.arange(-.4, 2, 0.1),
                constit_to_save = ['M2','S2','K2','N2','K1','O1','P1','Q1'], 
                semidiurnal_constit = ['M2','S2','K2','N2'],
                diurnal_constit = ['K1','O1','P1','Q1'],
                apply_ntr_filter = True ):
    '''
    Routine for analysis and comparison of model and observed SSH
    This routine calculates:
        1. Estimates of non-tidal residuals by subtracting an equivalent
           harmonic analysis, i.e. the same time period, missing data,
           constituent sets etc.
        2. Saves some selected harmonic information. Harmonic analysis is done
           using the utide package, with constituent sets based on the
           Rayleigh criterion
        3. NTR stats: MAE, RMSE, correlation, climatology
        4. SSH stats: climatology
        5. Tide stats: RMSE, MAE, correlation, climatology for semidiurnal
           and diurnal frequency bands.
        6. Threshold analysis of NTR. Sum of peaks, time, daily max, monthly
           max over specified thresholds.
    
    INPUTS
     fn_ext               : Filepath to extracted SSH from analyse_ssh
     fn_out               : Filepath to desired output analysis file
     thresholds           : Array of NTR thresholds for analysis
     constit_to_save      : List of constituents to save amplitudes/phase
     semidiurnal_constit  : List of constituents to include in semidiurnal band
     diurnal_constit      : List of constituents to include in diurnal band
     apply_ntr_filter     : If true, apply Savgol filter to non-tidal residuals
                            before analysis.
    '''
    min_datapoints=700 
    ds_ssh = xr.open_dataset(fn_ext) 
    
    # Define Dimension Sizes
    n_port = ds_ssh.dims['port']
    n_time = ds_ssh.dims['time']
    n_constit = len(constit_to_save)
    n_thresholds = len(thresholds)
    seasons = ['DJF','JJA','MAM','SON','All']
    n_seasons = len(seasons)
    freq_bands = ['diurnal', 'semidiurnal', 'all']
    n_freq_bands = len(freq_bands)
    
    # Remove flagged locations
    ds_ssh.ssh_mod[ds_ssh.bad_flag.values] = np.nan
    ds_ssh.ssh_obs[ds_ssh.bad_flag.values] = np.nan 
    
    # NTR dataset
    ds_ntr = xr.Dataset(coords = dict(
                            time = ('time', ds_ssh.time.values),
                            longitude = ('port', ds_ssh.longitude.values),
                            latitude = ('port', ds_ssh.latitude.values)),
                        data_vars = dict(
                            ntr_mod = (['port','time'], np.zeros((n_port, n_time))*np.nan),
                            ntr_obs = (['port','time'], np.zeros((n_port, n_time))*np.nan),
                            ntr_err = (['port','time'], np.zeros((n_port, n_time))*np.nan),
                            ntr_square_err = (['port','time'], np.zeros((n_port, n_time))*np.nan),
                            ntr_abs_err = (['port','time'], np.zeros((n_port, n_time))*np.nan)))
    
    ds_tide = xr.Dataset(coords = dict(
                            time = ('time', ds_ssh.time.values),
                            longitude = ('port', ds_ssh.longitude.values),
                            latitude = ('port', ds_ssh.latitude.values),
                            freq_band = ('freq_band', freq_bands)),
                        data_vars = dict(
                            tide_mod = (['port','freq_band','time'], np.zeros((n_port, n_freq_bands, n_time))*np.nan),
                            tide_obs = (['port','freq_band','time'], np.zeros((n_port, n_freq_bands, n_time))*np.nan),
                            tide_err = (['port','freq_band','time'], np.zeros((n_port, n_freq_bands, n_time))*np.nan),
                            tide_square_err = (['port','freq_band','time'], np.zeros((n_port, n_freq_bands, n_time))*np.nan),
                            tide_abs_err = (['port','freq_band','time'], np.zeros((n_port, n_freq_bands,n_time))*np.nan)))
    
    # ANALYSIS dataset
    ds_stats = xr.Dataset(coords = dict(
                    longitude = ('port', ds_ssh.longitude.values),
                    latitude = ('port', ds_ssh.latitude.values),
                    time = ('time', ds_ssh.time.values),
                    season = ('season', seasons),
                    constituent = ('constituent', constit_to_save),
                    threshold = ('threshold', thresholds)),
               data_vars = dict(
                    a_mod = (['port','constituent'], np.zeros((n_port, n_constit))*np.nan),
                    a_obs = (['port','constituent'], np.zeros((n_port, n_constit))*np.nan),
                    g_mod = (['port','constituent'], np.zeros((n_port, n_constit))*np.nan),
                    g_obs = (['port','constituent'], np.zeros((n_port, n_constit))*np.nan),
                    a_err = (['port','constituent'], np.zeros((n_port, n_constit))*np.nan),
                    g_err = (['port','constituent'], np.zeros((n_port, n_constit))*np.nan),
                    ssh_std_obs = (['port','season'], np.zeros((n_port, n_seasons))*np.nan),
                    ssh_std_mod = (['port','season'], np.zeros((n_port, n_seasons))*np.nan),
                    ssh_std_err = (['port','season'], np.zeros((n_port, n_seasons))*np.nan),
                    ntr_corr = (['port','season'],  np.zeros((n_port, n_seasons))*np.nan),
                    ntr_mae  = (['port','season'], np.zeros((n_port, n_seasons))*np.nan),
                    ntr_me  = (['port','season'], np.zeros((n_port, n_seasons))*np.nan),
                    ntr_rmse = (['port','season'],  np.zeros((n_port, n_seasons))*np.nan),
                    ntr_err_std = (['port','season'],  np.zeros((n_port, n_seasons))*np.nan),
                    ntr_std_obs = (['port','season'],   np.zeros((n_port, n_seasons))*np.nan),
                    ntr_std_mod = (['port','season'],   np.zeros((n_port, n_seasons))*np.nan),
                    thresh_peak_mod  = (['port', 'threshold'], np.zeros((n_port, n_thresholds))),
                    thresh_peak_obs  = (['port', 'threshold'], np.zeros((n_port, n_thresholds))),
                    thresh_time_mod = (['port', 'threshold'], np.zeros((n_port, n_thresholds))),
                    thresh_time_obs = (['port', 'threshold'], np.zeros((n_port, n_thresholds))),
                    thresh_dailymax_mod = (['port', 'threshold'], np.zeros((n_port, n_thresholds))),
                    thresh_dailymax_obs = (['port', 'threshold'], np.zeros((n_port, n_thresholds))),
                    thresh_monthlymax_mod = (['port', 'threshold'], np.zeros((n_port, n_thresholds))),
                    thresh_monthlymax_obs = (['port', 'threshold'], np.zeros((n_port, n_thresholds))))) 
    
    # Identify seasons
    month_season_dict = {1:0, 2:0, 3:2, 4:2, 5:2, 6:1,
                         7:1, 8:1, 9:3, 10:3, 11:3, 12:0}
    time = ds_ssh.time.values
    pd_time = pd.to_datetime(time)
    pd_month = pd_time.month
    pd_season = [month_season_dict[ii] for ii in pd_month]
    
    # Loop over tide gauge locations, perform analysis per location
    for pp in range(0,n_port):
        
        # Temporary in-loop datasets
        ds_ssh_port = ds_ssh.isel(port=pp).load()
        ssh_mod = ds_ssh_port.ssh_mod
        ssh_obs = ds_ssh_port.ssh_obs
        mask = ds_ssh.mask.isel(port=pp).load().values
        
        if all(np.isnan(ssh_mod)) or all(np.isnan(ssh_obs)):
            print('reject 1')
            continue
        
        if sum(~np.isnan(ssh_obs.values))<min_datapoints:
            print('reject 2')
            print(sum(~np.isnan(ssh_obs.values)))

        # Harmonic analysis datenums
        ha_time  = mdates.date2num(time)
        
        # Do harmonic analysis using UTide
        uts_obs = ut.solve(ha_time, ssh_obs.values, lat=ssh_obs.latitude.values)
        uts_mod = ut.solve(ha_time, ssh_mod.values, lat=ssh_mod.latitude.values)
        
        # Reconstruct full tidal signal 
        tide_obs = np.array( ut.reconstruct(ha_time, uts_obs).h)
        tide_mod = np.array( ut.reconstruct(ha_time, uts_mod).h)
        tide_obs[mask] = np.nan
        tide_mod[mask] = np.nan
        ds_tide['tide_obs'][pp, -1, :] = tide_obs
        ds_tide['tide_mod'][pp, -1, :] = tide_mod
        
        # Reconstruct partial semidiurnal tidal signal 
        tide_2_obs = np.array( ut.reconstruct(ha_time, uts_obs,
                                              constit = semidiurnal_constit).h)
        tide_2_mod = np.array( ut.reconstruct(ha_time, uts_mod, 
                                              constit = semidiurnal_constit).h)
        tide_2_obs[mask] = np.nan
        tide_2_mod[mask] = np.nan
        ds_tide['tide_obs'][pp, 1, :] = tide_2_obs
        ds_tide['tide_mod'][pp, 1, :] = tide_2_mod
        
        # # Reconstruct partial diurnal tidal signal 
        tide_1_obs = np.array( ut.reconstruct(ha_time, uts_obs,
                                              constit = diurnal_constit).h)
        tide_1_mod = np.array( ut.reconstruct(ha_time, uts_mod,
                                              constit = diurnal_constit).h)
        tide_1_obs[mask] = np.nan
        tide_1_mod[mask] = np.nan
        ds_tide['tide_obs'][pp, 0, :] = tide_1_obs
        ds_tide['tide_mod'][pp, 0, :] = tide_1_mod
        
        # TWL: SAVE constituents
        a_dict_obs = dict( zip(uts_obs.name, uts_obs.A) )
        a_dict_mod = dict( zip(uts_mod.name, uts_mod.A) )
        g_dict_obs = dict( zip(uts_obs.name, uts_obs.g) )
        g_dict_mod = dict( zip(uts_mod.name, uts_mod.g) )
        
        for cc in range(0, len(constit_to_save)):
            if constit_to_save[cc] in uts_obs.name:
                ds_stats['a_mod'][pp,cc] = a_dict_mod[constit_to_save[cc]] 
                ds_stats['a_obs'][pp,cc] = a_dict_obs[constit_to_save[cc]] 
                ds_stats['a_err'][pp,cc] = a_dict_mod[constit_to_save[cc]] - a_dict_obs[constit_to_save[cc]] 
                ds_stats['g_mod'][pp,cc] = g_dict_mod[constit_to_save[cc]] 
                ds_stats['g_obs'][pp,cc] = g_dict_obs[constit_to_save[cc]]
                ds_stats['g_err'][pp,cc] = g_dict_mod[constit_to_save[cc]] - g_dict_obs[constit_to_save[cc]] 
        
        # NTR: Calculate non tidal residuals
        ntr_obs = ssh_obs.values - tide_obs
        ntr_mod = ssh_mod.values - tide_mod
        
        # NTR: Apply filter if wanted
        if apply_ntr_filter:
            ntr_obs = signal.savgol_filter(ntr_obs,25,3)
            ntr_mod = signal.savgol_filter(ntr_mod,25,3)
            
        if sum(~np.isnan(ntr_obs)) < min_datapoints:
            continue
            
        ntr_err = ntr_mod - ntr_obs
        ds_ntr['ntr_obs'][pp] = ntr_obs
        ds_ntr['ntr_mod'][pp] = ntr_mod
        ds_ntr['ntr_err'][pp] = ntr_err
        ds_ntr['ntr_abs_err'][pp] = np.abs(ntr_err)
        ds_ntr['ntr_square_err'][pp] = ntr_err**2
        
        # Make masked arrays for seasonal correlation calculation
        ntr_obs = np.ma.masked_invalid(ntr_obs)
        ntr_mod = np.ma.masked_invalid(ntr_mod)
        ds_stats['ntr_corr'][pp,4] = np.ma.corrcoef(ntr_obs, ntr_mod)[1,0]
        for ss in range(0,4):
            season_ind = pd_season == ss
            if np.sum(season_ind)>100:
                tmp_obs = ntr_obs[season_ind]
                tmp_mod = ntr_mod[season_ind]
                ds_stats['ntr_corr'][pp,4] = np.ma.corrcoef(tmp_obs, tmp_mod)[1,0]
            
        # Identify NTR peaks for threshold analysis
        pk_ind_ntr_obs,_ = signal.find_peaks(ntr_obs, distance = 12)
        pk_ind_ntr_mod,_ = signal.find_peaks(ntr_mod, distance = 12)
        pk_ntr_obs = ntr_obs[pk_ind_ntr_obs]
        pk_ntr_mod = ntr_mod[pk_ind_ntr_mod]
        
        # Calculate daily and monthly maxima for threshold analysis
        ds_daily = ds_ntr.groupby('time.day')
        ds_daily_max = ds_daily.max(skipna=True)
        ds_monthly = ds_ntr.groupby('time.month')
        ds_monthly_max = ds_monthly.max(skipna=True)
        
        # Threshold Analysis
        for nn in range(0,n_thresholds):
            threshn = thresholds[nn]
            # NTR: Threshold Frequency (Peaks)
            ds_stats['thresh_peak_mod'][pp, nn] = np.sum( pk_ntr_mod >= threshn)
            ds_stats['thresh_peak_obs'][pp, nn] = np.sum( pk_ntr_obs >= threshn)
            
            # NTR: Threshold integral (Time over threshold)
            ds_stats['thresh_time_mod'][pp, nn] = np.sum( ntr_mod >= threshn)
            ds_stats['thresh_time_obs'][pp, nn] = np.sum( ntr_obs >=threshn)
            
            # NTR: Number of daily maxima over threshold
            ds_stats['thresh_dailymax_mod'][pp, nn] = np.sum( ds_daily_max.ntr_mod.values >= threshn)
            ds_stats['thresh_dailymax_obs'][pp, nn] = np.sum( ds_daily_max.ntr_obs.values >= threshn)
            
            # NTR: Number of monthly maxima over threshold
            ds_stats['thresh_monthlymax_mod'][pp, nn] = np.sum( ds_monthly_max.ntr_mod.values >= threshn)
            ds_stats['thresh_monthlymax_obs'][pp, nn] = np.sum( ds_monthly_max.ntr_mod.values >= threshn)
            
    
    # Seasonal Climatology
    ntr_seasonal = ds_ntr.groupby('time.season')
    ntr_seasonal_std = ntr_seasonal.std(skipna=True)
    ntr_seasonal_mean = ntr_seasonal.mean(skipna=True)
    ssh_seasonal = ds_ssh.groupby('time.season')
    ssh_seasonal_std = ssh_seasonal.std(skipna=True)
    
    sii = 0
    for ss in ntr_seasonal_std['season'].values:
        ind = seasons.index(ss)
        
        ds_stats['ntr_std_mod'][:, ind] = ntr_seasonal_std.ntr_mod.sel(season=ss)
        ds_stats['ntr_std_obs'][:, ind] = ntr_seasonal_std.ntr_obs.sel(season=ss)
        ds_stats['ntr_err_std'][:, ind] = ntr_seasonal_std.ntr_err.sel(season=ss)
        
        ds_stats['ntr_mae'][:, ind] = ntr_seasonal_mean.ntr_abs_err.sel(season=ss)
        ds_stats['ntr_rmse'][:, ind] = np.nanmean( ntr_seasonal_mean.ntr_square_err.sel(season=ss) )
        ds_stats['ntr_me'][:, ind] = ntr_seasonal_mean.ntr_err.sel(season=ss)
        
        ds_stats['ssh_std_mod'][:, ind] = ssh_seasonal_std.ssh_mod.sel(season=ss)
        ds_stats['ssh_std_obs'][:, ind] = ssh_seasonal_std.ssh_obs.sel(season=ss)
        ds_stats['ssh_std_err'][:, ind] = ssh_seasonal_std.ssh_mod.sel(season=ss) - ssh_seasonal_std.ssh_obs.sel(season=ss)
        sii+=1
        
    # Annual means and standard deviations
    ntr_std = ds_ntr.std(dim='time', skipna=True)
    ssh_std = ds_ssh.std(dim='time', skipna=True)
    ntr_mean = ds_ntr.mean(dim='time', skipna=True)
    
    ds_stats['ntr_std_mod'][:, 4] = ntr_std.ntr_mod
    ds_stats['ntr_std_obs'][:, 4] = ntr_std.ntr_obs
    ds_stats['ntr_err_std'][:, 4] = ntr_std.ntr_err
    
    ds_stats['ntr_mae'][:, 4] = ntr_mean.ntr_abs_err
    ds_stats['ntr_rmse'][:, 4] = np.nanmean( ntr_mean.ntr_square_err )
    ds_stats['ntr_me'][:, 4] = ntr_mean.ntr_err
    
    ds_stats['ssh_std_mod'][:, 4] = ssh_std.ssh_mod
    ds_stats['ssh_std_obs'][:, 4] = ssh_std.ssh_obs
    ds_stats['ssh_std_err'][:, 4] = ssh_std.ssh_mod - ssh_std.ssh_obs
    
    ds_stats = xr.merge((ds_ssh, ds_ntr, ds_stats, ds_tide))
    
    write_ds_to_file(ds_stats, fn_out)
    
def analyse_ssh(fn_ext, fn_out, thresholds = np.arange(-.4, 2, 0.1),
                constit_to_save = ['M2','S2','K2','N2','K1','O1','P1','Q1'], 
                semidiurnal_constit = ['M2','S2','K2','N2'],
                diurnal_constit = ['K1','O1','P1','Q1'],
                apply_ntr_filter = True ):
    '''
    Routine for analysis and comparison of model and observed SSH
    This routine calculates:
        1. Estimates of non-tidal residuals by subtracting an equivalent
           harmonic analysis, i.e. the same time period, missing data,
           constituent sets etc.
        2. Saves some selected harmonic information. Harmonic analysis is done
           using the utide package, with constituent sets based on the
           Rayleigh criterion
        3. NTR stats: MAE, RMSE, correlation, climatology
        4. SSH stats: climatology
        5. Tide stats: RMSE, MAE, correlation, climatology for semidiurnal
           and diurnal frequency bands.
        6. Threshold analysis of NTR. Sum of peaks, time, daily max, monthly
           max over specified thresholds.
    
    INPUTS
     fn_ext               : Filepath to extracted SSH from analyse_ssh
     fn_out               : Filepath to desired output analysis file
     thresholds           : Array of NTR thresholds for analysis
     constit_to_save      : List of constituents to save amplitudes/phase
     semidiurnal_constit  : List of constituents to include in semidiurnal band
     diurnal_constit      : List of constituents to include in diurnal band
     apply_ntr_filter     : If true, apply Savgol filter to non-tidal residuals
                            before analysis.
    '''
    min_datapoints=700 
    ds_ssh = xr.open_dataset(fn_ext) 
    
    # Define Dimension Sizes
    n_port = ds_ssh.dims['port']
    n_time = ds_ssh.dims['time']
    n_constit = len(constit_to_save)
    n_thresholds = len(thresholds)
    seasons = ['DJF','JJA','MAM','SON','All']
    n_seasons = len(seasons)
    freq_bands = ['diurnal', 'semidiurnal', 'all']
    n_freq_bands = len(freq_bands)
    
    # Remove flagged locations
    ds_ssh.ssh_mod[ds_ssh.bad_flag.values] = np.nan
    ds_ssh.ssh_obs[ds_ssh.bad_flag.values] = np.nan 
    
    # NTR dataset
    ds_ntr = xr.Dataset(coords = dict(
                            time = ('time', ds_ssh.time.values),
                            longitude = ('port', ds_ssh.longitude.values),
                            latitude = ('port', ds_ssh.latitude.values)),
                        data_vars = dict(
                            ntr_mod = (['port','time'], np.zeros((n_port, n_time))*np.nan),
                            ntr_obs = (['port','time'], np.zeros((n_port, n_time))*np.nan),
                            ntr_err = (['port','time'], np.zeros((n_port, n_time))*np.nan),
                            ntr_square_err = (['port','time'], np.zeros((n_port, n_time))*np.nan),
                            ntr_abs_err = (['port','time'], np.zeros((n_port, n_time))*np.nan)))
    
    ds_tide = xr.Dataset(coords = dict(
                            time = ('time', ds_ssh.time.values),
                            longitude = ('port', ds_ssh.longitude.values),
                            latitude = ('port', ds_ssh.latitude.values),
                            freq_band = ('freq_band', freq_bands)),
                        data_vars = dict(
                            tide_mod = (['port','freq_band','time'], np.zeros((n_port, n_freq_bands, n_time))*np.nan),
                            tide_obs = (['port','freq_band','time'], np.zeros((n_port, n_freq_bands, n_time))*np.nan),
                            tide_err = (['port','freq_band','time'], np.zeros((n_port, n_freq_bands, n_time))*np.nan),
                            tide_square_err = (['port','freq_band','time'], np.zeros((n_port, n_freq_bands, n_time))*np.nan),
                            tide_abs_err = (['port','freq_band','time'], np.zeros((n_port, n_freq_bands,n_time))*np.nan)))
    
    # ANALYSIS dataset
    ds_stats = xr.Dataset(coords = dict(
                    longitude = ('port', ds_ssh.longitude.values),
                    latitude = ('port', ds_ssh.latitude.values),
                    time = ('time', ds_ssh.time.values),
                    season = ('season', seasons),
                    constituent = ('constituent', constit_to_save),
                    threshold = ('threshold', thresholds)),
               data_vars = dict(
                    a_mod = (['port','constituent'], np.zeros((n_port, n_constit))*np.nan),
                    a_obs = (['port','constituent'], np.zeros((n_port, n_constit))*np.nan),
                    g_mod = (['port','constituent'], np.zeros((n_port, n_constit))*np.nan),
                    g_obs = (['port','constituent'], np.zeros((n_port, n_constit))*np.nan),
                    a_err = (['port','constituent'], np.zeros((n_port, n_constit))*np.nan),
                    g_err = (['port','constituent'], np.zeros((n_port, n_constit))*np.nan),
                    ssh_std_obs = (['port','season'], np.zeros((n_port, n_seasons))*np.nan),
                    ssh_std_mod = (['port','season'], np.zeros((n_port, n_seasons))*np.nan),
                    ssh_std_err = (['port','season'], np.zeros((n_port, n_seasons))*np.nan),
                    ntr_corr = (['port','season'],  np.zeros((n_port, n_seasons))*np.nan),
                    ntr_mae  = (['port','season'], np.zeros((n_port, n_seasons))*np.nan),
                    ntr_me  = (['port','season'], np.zeros((n_port, n_seasons))*np.nan),
                    ntr_rmse = (['port','season'],  np.zeros((n_port, n_seasons))*np.nan),
                    ntr_err_std = (['port','season'],  np.zeros((n_port, n_seasons))*np.nan),
                    ntr_std_obs = (['port','season'],   np.zeros((n_port, n_seasons))*np.nan),
                    ntr_std_mod = (['port','season'],   np.zeros((n_port, n_seasons))*np.nan),
                    thresh_peak_mod  = (['port', 'threshold'], np.zeros((n_port, n_thresholds))),
                    thresh_peak_obs  = (['port', 'threshold'], np.zeros((n_port, n_thresholds))),
                    thresh_time_mod = (['port', 'threshold'], np.zeros((n_port, n_thresholds))),
                    thresh_time_obs = (['port', 'threshold'], np.zeros((n_port, n_thresholds))),
                    thresh_dailymax_mod = (['port', 'threshold'], np.zeros((n_port, n_thresholds))),
                    thresh_dailymax_obs = (['port', 'threshold'], np.zeros((n_port, n_thresholds))),
                    thresh_monthlymax_mod = (['port', 'threshold'], np.zeros((n_port, n_thresholds))),
                    thresh_monthlymax_obs = (['port', 'threshold'], np.zeros((n_port, n_thresholds))))) 
    
    # Identify seasons
    month_season_dict = {1:0, 2:0, 3:2, 4:2, 5:2, 6:1,
                         7:1, 8:1, 9:3, 10:3, 11:3, 12:0}
    time = ds_ssh.time.values
    pd_time = pd.to_datetime(time)
    pd_month = pd_time.month
    pd_season = [month_season_dict[ii] for ii in pd_month]
    
    # Loop over tide gauge locations, perform analysis per location
    for pp in range(0,n_port):
        
        # Temporary in-loop datasets
        ds_ssh_port = ds_ssh.isel(port=pp).load()
        ssh_mod = ds_ssh_port.ssh_mod
        ssh_obs = ds_ssh_port.ssh_obs
        mask = ds_ssh.mask.isel(port=pp).load().values
        
        if all(np.isnan(ssh_mod)) or all(np.isnan(ssh_obs)):
            print('reject 1')
            continue
        
        if sum(~np.isnan(ssh_obs.values))<min_datapoints:
            print('reject 2')
            print(sum(~np.isnan(ssh_obs.values)))

        # Harmonic analysis datenums
        ha_time  = mdates.date2num(time)
        
        # Do harmonic analysis using UTide
        uts_obs = ut.solve(ha_time, ssh_obs.values, lat=ssh_obs.latitude.values)
        uts_mod = ut.solve(ha_time, ssh_mod.values, lat=ssh_mod.latitude.values)
        
        # Reconstruct full tidal signal 
        tide_obs = np.array( ut.reconstruct(ha_time, uts_obs).h)
        tide_mod = np.array( ut.reconstruct(ha_time, uts_mod).h)
        tide_obs[mask] = np.nan
        tide_mod[mask] = np.nan
        ds_tide['tide_obs'][pp, -1, :] = tide_obs
        ds_tide['tide_mod'][pp, -1, :] = tide_mod
        
        # Reconstruct partial semidiurnal tidal signal 
        tide_2_obs = np.array( ut.reconstruct(ha_time, uts_obs,
                                              constit = semidiurnal_constit).h)
        tide_2_mod = np.array( ut.reconstruct(ha_time, uts_mod, 
                                              constit = semidiurnal_constit).h)
        tide_2_obs[mask] = np.nan
        tide_2_mod[mask] = np.nan
        ds_tide['tide_obs'][pp, 1, :] = tide_2_obs
        ds_tide['tide_mod'][pp, 1, :] = tide_2_mod
        
        # # Reconstruct partial diurnal tidal signal 
        tide_1_obs = np.array( ut.reconstruct(ha_time, uts_obs,
                                              constit = diurnal_constit).h)
        tide_1_mod = np.array( ut.reconstruct(ha_time, uts_mod,
                                              constit = diurnal_constit).h)
        tide_1_obs[mask] = np.nan
        tide_1_mod[mask] = np.nan
        ds_tide['tide_obs'][pp, 0, :] = tide_1_obs
        ds_tide['tide_mod'][pp, 0, :] = tide_1_mod
        
        # TWL: SAVE constituents
        a_dict_obs = dict( zip(uts_obs.name, uts_obs.A) )
        a_dict_mod = dict( zip(uts_mod.name, uts_mod.A) )
        g_dict_obs = dict( zip(uts_obs.name, uts_obs.g) )
        g_dict_mod = dict( zip(uts_mod.name, uts_mod.g) )
        
        for cc in range(0, len(constit_to_save)):
            if constit_to_save[cc] in uts_obs.name:
                ds_stats['a_mod'][pp,cc] = a_dict_mod[constit_to_save[cc]] 
                ds_stats['a_obs'][pp,cc] = a_dict_obs[constit_to_save[cc]] 
                ds_stats['a_err'][pp,cc] = a_dict_mod[constit_to_save[cc]] - a_dict_obs[constit_to_save[cc]] 
                ds_stats['g_mod'][pp,cc] = g_dict_mod[constit_to_save[cc]] 
                ds_stats['g_obs'][pp,cc] = g_dict_obs[constit_to_save[cc]]
                ds_stats['g_err'][pp,cc] = g_dict_mod[constit_to_save[cc]] - g_dict_obs[constit_to_save[cc]] 
        
        # NTR: Calculate non tidal residuals
        ntr_obs = ssh_obs.values - tide_obs
        ntr_mod = ssh_mod.values - tide_mod
        
        # NTR: Apply filter if wanted
        if apply_ntr_filter:
            ntr_obs = signal.savgol_filter(ntr_obs,25,3)
            ntr_mod = signal.savgol_filter(ntr_mod,25,3)
            
        if sum(~np.isnan(ntr_obs)) < min_datapoints:
            continue
            
        ntr_err = ntr_mod - ntr_obs
        ds_ntr['ntr_obs'][pp] = ntr_obs
        ds_ntr['ntr_mod'][pp] = ntr_mod
        ds_ntr['ntr_err'][pp] = ntr_err
        ds_ntr['ntr_abs_err'][pp] = np.abs(ntr_err)
        ds_ntr['ntr_square_err'][pp] = ntr_err**2
        
        # Make masked arrays for seasonal correlation calculation
        ntr_obs = np.ma.masked_invalid(ntr_obs)
        ntr_mod = np.ma.masked_invalid(ntr_mod)
        ds_stats['ntr_corr'][pp,4] = np.ma.corrcoef(ntr_obs, ntr_mod)[1,0]
        for ss in range(0,4):
            season_ind = pd_season == ss
            if np.sum(season_ind)>100:
                tmp_obs = ntr_obs[season_ind]
                tmp_mod = ntr_mod[season_ind]
                ds_stats['ntr_corr'][pp,4] = np.ma.corrcoef(tmp_obs, tmp_mod)[1,0]
            
        # Identify NTR peaks for threshold analysis
        pk_ind_ntr_obs,_ = signal.find_peaks(ntr_obs, distance = 12)
        pk_ind_ntr_mod,_ = signal.find_peaks(ntr_mod, distance = 12)
        pk_ntr_obs = ntr_obs[pk_ind_ntr_obs]
        pk_ntr_mod = ntr_mod[pk_ind_ntr_mod]
        
        # Calculate daily and monthly maxima for threshold analysis
        ds_daily = ds_ntr.groupby('time.day')
        ds_daily_max = ds_daily.max(skipna=True)
        ds_monthly = ds_ntr.groupby('time.month')
        ds_monthly_max = ds_monthly.max(skipna=True)
        
        # Threshold Analysis
        for nn in range(0,n_thresholds):
            threshn = thresholds[nn]
            # NTR: Threshold Frequency (Peaks)
            ds_stats['thresh_peak_mod'][pp, nn] = np.sum( pk_ntr_mod >= threshn)
            ds_stats['thresh_peak_obs'][pp, nn] = np.sum( pk_ntr_obs >= threshn)
            
            # NTR: Threshold integral (Time over threshold)
            ds_stats['thresh_time_mod'][pp, nn] = np.sum( ntr_mod >= threshn)
            ds_stats['thresh_time_obs'][pp, nn] = np.sum( ntr_obs >=threshn)
            
            # NTR: Number of daily maxima over threshold
            ds_stats['thresh_dailymax_mod'][pp, nn] = np.sum( ds_daily_max.ntr_mod.values >= threshn)
            ds_stats['thresh_dailymax_obs'][pp, nn] = np.sum( ds_daily_max.ntr_obs.values >= threshn)
            
            # NTR: Number of monthly maxima over threshold
            ds_stats['thresh_monthlymax_mod'][pp, nn] = np.sum( ds_monthly_max.ntr_mod.values >= threshn)
            ds_stats['thresh_monthlymax_obs'][pp, nn] = np.sum( ds_monthly_max.ntr_mod.values >= threshn)
            
    
    # Seasonal Climatology
    ntr_seasonal = ds_ntr.groupby('time.season')
    ntr_seasonal_std = ntr_seasonal.std(skipna=True)
    ntr_seasonal_mean = ntr_seasonal.mean(skipna=True)
    ssh_seasonal = ds_ssh.groupby('time.season')
    ssh_seasonal_std = ssh_seasonal.std(skipna=True)
    
    sii = 0
    for ss in ntr_seasonal_std['season'].values:
        ind = seasons.index(ss)
        
        ds_stats['ntr_std_mod'][:, ind] = ntr_seasonal_std.ntr_mod.sel(season=ss)
        ds_stats['ntr_std_obs'][:, ind] = ntr_seasonal_std.ntr_obs.sel(season=ss)
        ds_stats['ntr_err_std'][:, ind] = ntr_seasonal_std.ntr_err.sel(season=ss)
        
        ds_stats['ntr_mae'][:, ind] = ntr_seasonal_mean.ntr_abs_err.sel(season=ss)
        ds_stats['ntr_rmse'][:, ind] = np.nanmean( ntr_seasonal_mean.ntr_square_err.sel(season=ss) )
        ds_stats['ntr_me'][:, ind] = ntr_seasonal_mean.ntr_err.sel(season=ss)
        
        ds_stats['ssh_std_mod'][:, ind] = ssh_seasonal_std.ssh_mod.sel(season=ss)
        ds_stats['ssh_std_obs'][:, ind] = ssh_seasonal_std.ssh_obs.sel(season=ss)
        ds_stats['ssh_std_err'][:, ind] = ssh_seasonal_std.ssh_mod.sel(season=ss) - ssh_seasonal_std.ssh_obs.sel(season=ss)
        sii+=1
        
    # Annual means and standard deviations
    ntr_std = ds_ntr.std(dim='time', skipna=True)
    ssh_std = ds_ssh.std(dim='time', skipna=True)
    ntr_mean = ds_ntr.mean(dim='time', skipna=True)
    
    ds_stats['ntr_std_mod'][:, 4] = ntr_std.ntr_mod
    ds_stats['ntr_std_obs'][:, 4] = ntr_std.ntr_obs
    ds_stats['ntr_err_std'][:, 4] = ntr_std.ntr_err
    
    ds_stats['ntr_mae'][:, 4] = ntr_mean.ntr_abs_err
    ds_stats['ntr_rmse'][:, 4] = np.nanmean( ntr_mean.ntr_square_err )
    ds_stats['ntr_me'][:, 4] = ntr_mean.ntr_err
    
    ds_stats['ssh_std_mod'][:, 4] = ssh_std.ssh_mod
    ds_stats['ssh_std_obs'][:, 4] = ssh_std.ssh_obs
    ds_stats['ssh_std_err'][:, 4] = ssh_std.ssh_mod - ssh_std.ssh_obs
    
    ds_stats = xr.merge((ds_ssh, ds_ntr, ds_stats, ds_tide))
    
    write_ds_to_file(ds_stats, fn_out)
    
def extract_ssh(fn_nemo_data, fn_nemo_domain, fn_nemo_cfg, fn_obs, fn_out,
                     chunks = {'time_counter':100}, dist_omit = 5):
                     
    '''
    Routine for extraction of model ssh at obs locations.
    The tidegauge file should be a netcdf file with dimension
    ('port', 'time') and variables 'ssh', 'time', 'longitude', 'latitude'.
    All ports are expected to have data on a common frequency, so some
    preprocessing of obs is required.
    
    INPUTS
     fn_nemo_data    : Absoltue path to NEMO data file(s)
     fn_nemo_domain  : Absolute path to NEMO domain file
     fn_obs          : Absolute path to Tidegauge data file
     fn_out          : Absolute path to output file
     chunks          : xarray chunking dictionary
     dist_omit       : Distance (km) from obs point to reject nearest model point
    '''
    
    # Read NEMO data
    nemo = coast.Gridded(fn_nemo_data, fn_nemo_domain, config=fn_nemo_cfg,
                              multiple=True).dataset
    
    # Get NEMO landmask
    landmask = np.array(nemo.bottom_level.values.squeeze() == 0)
    
    # Read OBS data from tg file
    obs = xr.open_dataset(fn_obs)
    
    # Subset obs by NEMO domain
    lonmax = np.nanmax(nemo.longitude)
    lonmin = np.nanmin(nemo.longitude)
    latmax = np.nanmax(nemo.latitude)
    latmin = np.nanmin(nemo.latitude)
    ind = gu.subset_indices_lonlat_box(obs.longitude, obs.latitude, 
                                       lonmin, lonmax, latmin, latmax)
    obs = obs.isel(port=ind[0])
    
    # Determine spatial indices
    ind2D = gu.nearest_indices_2d(nemo.longitude, nemo.latitude,
                                obs.longitude, obs.latitude,
                                mask = landmask)
    
    # Extract spatial time series
    nemo_extracted = nemo.isel(x_dim = ind2D[0], y_dim = ind2D[1])
    nemo_extracted = nemo_extracted.swap_dims({'dim_0':'port'})
    
    # Compute data (takes a while..)
    with ProgressBar():
        nemo_extracted.load()
    
    # Check interpolation distances
    n_port = nemo_extracted.dims['port']
    bad_flag = np.zeros(n_port, dtype=bool)
    interp_dist = gu.calculate_haversine_distance(nemo_extracted.longitude, 
                                                  nemo_extracted.latitude, 
                                                  obs.longitude.values,
                                                  obs.latitude.values)
    omit_ind = interp_dist>dist_omit
    bad_flag[omit_ind] = True
    
    # Align timings
    obs = obs.interp(time = nemo_extracted.time_instant.values, method = 'linear')
    
    # Apply shared mask
    ssh_mod = nemo_extracted.ssh.values.T
    ssh_obs = obs.ssh.values
    mask_mod = np.isnan(ssh_mod)
    mask_obs = np.isnan(ssh_obs)
    shared_mask = np.logical_or(mask_mod, mask_obs)
    
    ssh_mod[shared_mask] = np.nan
    ssh_obs[shared_mask] = np.nan
    
    ext = xr.Dataset(coords = dict(
                    longitude = ('port', obs.longitude.values),
                    latitude = ('port', obs.latitude.values),
                    time = ('time', nemo_extracted.time_instant.values)),
               data_vars = dict(
                    ssh_mod = (['port','time'], ssh_mod),
                    ssh_obs  = (['port','time'], ssh_obs),
                    bad_flag = (['port'], bad_flag),
                    mask = (['port','time'], shared_mask)))
    
    write_ds_to_file(ext, fn_out)
        
class plot_single_cfg():
    
    def __init__(self, fn_ssh_hourly_stats, dn_out, run_name, file_type='.png'):
    
        stats = xr.open_dataset(fn_ssh_hourly_stats, engine="netcdf4")
        
        print(f"stats:{stats}")
        #print(f"stats.amp_err:{stats.amp_err}")
        lonmax = np.max(stats.longitude)
        lonmin = np.min(stats.longitude)
        latmax = np.max(stats.latitude)
        latmin = np.min(stats.latitude)
        lonbounds = [lonmin-4, lonmax+4]
        latbounds = [latmin-4, latmax+4]
        
        ### GEOGRAPHICAL SCATTER PLOTS
        # Plot correlations
        f,a = pu.create_geo_axes(lonbounds, latbounds)
        sca = a.scatter(stats.longitude, stats.latitude, c=stats.ntr_corr.sel(season="All"), 
                        vmin=.85, vmax=1,
                  edgecolors='k', linewidths=.5, zorder=100)
        f.colorbar(sca)
        a.set_title('Hourly NTR Correlations with Tide Gauge | {0}'.format(run_name), fontsize=9)
        fn = "ssh_hourly_ntr_correlations_"+run_name+file_type
        f.savefig(os.path.join(dn_out, fn))
        plt.close('all')
        
        # Plot std_err
        f,a = pu.create_geo_axes(lonbounds, latbounds)
        sca = a.scatter(stats.longitude, stats.latitude, c=stats.ssh_std_err.sel(season="All"), 
                        vmin=-.15, vmax=.15,
                  edgecolors='k', linewidths=.5, zorder=100, cmap='seismic')
        f.colorbar(sca)
        a.set_title('Hourly SSH St. Dev Error (m) | {0}'.format(run_name), fontsize=9)
        fn = "ssh_hourly_stderr_"+run_name+file_type
        f.savefig(os.path.join(dn_out, fn))
        plt.close('all')
        
        # Plot mae
        f,a = pu.create_geo_axes(lonbounds, latbounds)
        sca = a.scatter(stats.longitude, stats.latitude, c=stats.ntr_mae.sel(season="All"), 
                        vmin=-.05, vmax=.05,
                  edgecolors='k', linewidths=.5, zorder=100)
        f.colorbar(sca)
        a.set_title('Hourly NTR MAE (m) | {0}'.format(run_name), fontsize=9)
        fn = "ssh_hourly_ntr_mae_"+run_name+file_type
        f.savefig(os.path.join(dn_out, fn))
        plt.close('all')
        
        ### HARMONIC ERRORS 
        constit = stats.constituent.values
        n_constit = len(constit)
        
        for cc in range(0,n_constit):
            # AMPLITUDE MAP
            f,a = pu.create_geo_axes(lonbounds, latbounds)
            sca = a.scatter(stats.longitude, stats.latitude, c=stats.a_err[:,cc], 
                            vmin=-.2, vmax=.2,
                      edgecolors='k', linewidths=.5, zorder=100, cmap='seismic')
            f.colorbar(sca)
            a.set_title('Amplitude Error (m) | {1} - Obs | {0}'.format(constit[cc], run_name), fontsize=9)
            fn = "ssh_hourly_amp_err_{0}_{1}{2}".format(constit[cc],run_name,file_type)
            f.savefig(os.path.join(dn_out, fn))
            plt.close('all')
            
            # PHASE MAP
            f,a = pu.create_geo_axes(lonbounds, latbounds)
            sca = a.scatter(stats.longitude, stats.latitude, c=stats.g_err[:,cc], 
                            vmin=-15, vmax=15,
                      edgecolors='k', linewidths=.5, zorder=100, cmap='seismic')
            f.colorbar(sca)
            a.set_title('Phase Error (deg) | {1} - Obs | {0}'.format(constit[cc], run_name), fontsize=9)
            fn = "ssh_hourly_pha_err_{0}_{1}{2}".format(constit[cc],run_name,file_type)
            f.savefig(os.path.join(dn_out, fn))
            plt.close('all')
            
            # AMPLITUDE SCATTER
            f,a = pu.scatter_with_fit(stats.a_mod[:,cc], stats.a_obs[:,cc])
            a.set_title('Amplitude Comparison | {0} | {1}'.format(constit[cc], run_name), fontsize=9)
            a.set_xlabel("Model Amplitude (m)")
            a.set_ylabel("Observed Amplitude (m)")
            fn = "ssh_hourly_scatter_amp_{0}_{1}{2}".format(constit[cc], run_name, file_type)
            f.savefig(os.path.join(dn_out, fn))
            plt.close('all')
            
            # PHASE SCATTER
            f,a = pu.scatter_with_fit(stats.g_mod[:,cc], stats.g_obs[:,cc])
            a.set_title('Phase Comparison | {0} | {1}'.format(constit[cc], run_name), fontsize=9)
            a.set_xlabel("Model Phase(deg)")
            a.set_ylabel("Observed Phase (deg)")
            fn = "ssh_hourly_scatter_pha_{0}_{1}{2}".format(constit[cc], run_name, file_type)
            f.savefig(os.path.join(dn_out, fn))
            plt.close('all')
        
        ### CLIM VAR ERROR
        # Jan
        m_ind = [0,3,6,9]
        m_str = ['Jan','Mar','Jul','Oct']
        n_m = len(m_ind)
        for mm in range(0,n_m):
            f,a = pu.create_geo_axes(lonbounds, latbounds)
            sca = a.scatter(stats.longitude, stats.latitude, c=stats.ntr_mod_clim_var[:,m_ind[mm]] - stats.ntr_obs_clim_var[:,m_ind[mm]], 
                            vmin=-.05, vmax=.05,
                      edgecolors='k', linewidths=.5, zorder=100, cmap='seismic')
            f.colorbar(sca)
            a.set_title('Climatological Variance Error | {1} | var({0}) - var(Obs)'.format(run_name, m_str[mm]), fontsize=9)
            fn = "ssh_hourly_clim_var_error_{0}_{1}{2}".format(m_str[mm],run_name,file_type)
            f.savefig(os.path.join(dn_out, fn))
            plt.close('all')
        
        ### FITTED SCATTER PLOTS
        # Plot st.dev comparisons
        f,a = pu.scatter_with_fit(stats.skew_mod.values.flatten(), 
                                  stats.skew_obs.values.flatten())
        a.set_title('Skew Surge Comparison | {0}'.format(run_name), fontsize=9)
        a.set_xlabel("Model Skew Surge (m)".format(run_name))
        a.set_ylabel("Observed Skew Surge (m)")
        a.set_xlim(-3,3)
        a.set_ylim(-3,3)
        a.set_title("Comparison of modelled and observed skew surges | {0}".format(run_name))
        fn = "ssh_hourly_skew_scatter_{0}{1}".format(run_name, file_type)
        f.savefig(os.path.join(dn_out, fn))
        plt.close('all')
        
        # THRESHOLD PLOTS
        prop = stats.thresh_freq_ntr_mod/stats.thresh_freq_ntr_obs
        plot_thresh = [5, 10, 15]
        for pp in range(0, len(plot_thresh)):
            prop_tmp = prop.isel(threshold=plot_thresh[pp])
            f,a = pu.create_geo_axes(lonbounds, latbounds)
            sca = a.scatter(stats.longitude, stats.latitude, c=prop_tmp, 
                            vmin=0, vmax=2,
                      edgecolors='k', linewidths=.5, zorder=100, cmap='seismic')
            f.colorbar(sca)
            a.set_title('Number of model NTR peaks as a proportion of observed NTR peaks \n Threshold = {0}m | {1}'.format(prop_tmp.threshold.values, run_name), fontsize=9)
            fn = "thresh_freq_{0}_{1}{2}".format(prop_tmp.threshold.values, run_name, file_type)
            f.savefig(os.path.join(dn_out, fn))
            plt.close('all')
            
        prop = stats.thresh_int_ntr_mod/stats.thresh_int_ntr_obs
        plot_thresh = [5, 10, 15]
        for pp in range(0, len(plot_thresh)):
            prop_tmp = prop.isel(threshold=plot_thresh[pp])
            f,a = pu.create_geo_axes(lonbounds, latbounds)
            sca = a.scatter(stats.longitude, stats.latitude, c=prop_tmp, 
                            vmin=0, vmax=2,
                      edgecolors='k', linewidths=.5, zorder=100, cmap='seismic')
            f.colorbar(sca)
            a.set_title('Model hours over threshold as a proportion of observed time \n Threshold = {0}m | {1}'.format(prop_tmp.threshold.values, run_name), fontsize=9)
            fn = "thresh_int_{0}_{1}{2}".format(prop_tmp.threshold.values, run_name, file_type)
            f.savefig(os.path.join(dn_out, fn))
            plt.close('all')
        
class plot_stats_ssh_hourly_multi_cfg():
    def __init__(self, fn_ssh_hourly_stats, dn_out, run_name, file_type='.png'):
        
        if type(fn_ssh_hourly_stats) is not list:
            fn_ssh_hourly_stats = [fn_ssh_hourly_stats]
            run_name = [run_name]
        
        if len(run_name) != len(fn_ssh_hourly_stats):
            print('Error: run_name length must be equal to fn_ssh_hourly_stats length.')
            return
        
        n_cfgs = len(fn_ssh_hourly_stats)

        stats = [xr.open_dataset(fn) for fn in fn_ssh_hourly_stats]

        ### THRESHOLD ANALYSIS
        
        # Ports Sum - Proportional
        sum_mod = [ss.thresh_freq_ntr_mod.sum(dim='port') for ss in stats]
        sum_obs = [ss.thresh_freq_ntr_obs.sum(dim='port') for ss in stats]
        xmax = np.nanmax(stats[0].threshold)
        f = plt.figure()
        for ss in range(0,n_cfgs):
            plt.plot(stats[ss].threshold, sum_mod[ss]/sum_obs[ss])
        plt.legend(run_name)
        plt.plot([-1000, 1000], [1,1], c='k', linestyle='--', linewidth=.75)
        plt.xlabel('Threshold (m)')
        plt.ylabel('Number of NTR peaks over threshold')
        plt.ylim(0,2)
        plt.xlim(0, xmax)
        plt.grid()
        plt.title('Proportional count of NTR peaks above threshold (Model/Obs)', fontsize=9)
        fn = "ntr_thresh_freq_prop"+file_type
        f.savefig(os.path.join(dn_out, fn))
        plt.close('all')
        
        # Ports Sum - Proportional INTEGRAL
        sum_mod = [ss.thresh_int_ntr_mod.sum(dim='port') for ss in stats]
        sum_obs = [ss.thresh_int_ntr_obs.sum(dim='port') for ss in stats]
        f = plt.figure()
        for ss in range(0,n_cfgs):
            plt.plot(stats[ss].threshold, sum_mod[ss]/sum_obs[ss])
        plt.legend(run_name)
        plt.plot([-1000, 1000],[1,1], c='k', linestyle='--', linewidth=.75)
        plt.xlabel('Threshold (m)')
        plt.ylabel('Time over threshold (hours)')
        plt.ylim(0,2)
        plt.xlim(0, xmax)
        plt.grid()
        plt.title('Proportional NTR integral over threshold (Model/Obs)', fontsize=9)
        fn = "ntr_thresh_int_prop"+file_type
        f.savefig(os.path.join(dn_out, fn))
        plt.close('all')
        
        
class plot_stats_ssh_hourly_compare_cfgs():
    def __init__(self, fn_ssh_hourly_stats1, fn_ssh_hourly_stats2, dn_out, 
                 run_name1, run_name2, file_type='.png'):
    
        stats1 = xr.open_dataset(fn_ssh_hourly_stats1)
        stats2 = xr.open_dataset(fn_ssh_hourly_stats2)
        
        lonmax = np.max(stats1.longitude)
        lonmin = np.min(stats1.longitude)
        latmax = np.max(stats1.latitude)
        latmin = np.min(stats1.latitude)
        lonbounds = [lonmin-4, lonmax+4]
        latbounds = [latmin-4, latmax+4]
        
        ### GEOGRAPHICAL SCATTER PLOTS
        # Plot correlations
        f,a = pu.create_geo_axes(lonbounds, latbounds)
        sca = a.scatter(stats.longitude, stats.latitude, c=stats.ntr_corr, 
                        vmin=.85, vmax=1,
                  edgecolors='k', linewidths=.5, zorder=100)
        f.colorbar(sca)
        a.set_title('Hourly NTR Correlations with Tide Gauge | {0}'.format(run_name), fontsize=9)
        fn = "ssh_hourly_ntr_correlations_"+run_name+file_type
        f.savefig(os.path.join(dn_out, fn))
        plt.close('all')
        
        # Plot std_err
        f,a = pu.create_geo_axes(lonbounds, latbounds)
        sca = a.scatter(stats.longitude, stats.latitude, c=stats.std_err, 
                        vmin=-.15, vmax=.15,
                  edgecolors='k', linewidths=.5, zorder=100, cmap='seismic')
        f.colorbar(sca)
        a.set_title('Hourly SSH St. Dev Error (m) | {0}'.format(run_name), fontsize=9)
        fn = "ssh_hourly_stderr_"+run_name+file_type
        f.savefig(os.path.join(dn_out, fn))
        plt.close('all')
        
        # Plot mae
        f,a = pu.create_geo_axes(lonbounds, latbounds)
        sca = a.scatter(stats.longitude, stats.latitude, c=stats.ntr_mae, 
                        vmin=-.05, vmax=.05,
                  edgecolors='k', linewidths=.5, zorder=100)
        f.colorbar(sca)
        a.set_title('Hourly NTR MAE (m) | {0}'.format(run_name), fontsize=9)
        fn = "ssh_hourly_ntr_mae_"+run_name+file_type
        f.savefig(os.path.join(dn_out, fn))
        plt.close('all')
        
        ### HARMONIC ERRORS 
        constit = stats.constituent.values
        n_constit = len(constit)
        
        for cc in range(0,n_constit):
            # AMPLITUDE MAP
            f,a = pu.create_geo_axes(lonbounds, latbounds)
            sca = a.scatter(stats.longitude, stats.latitude, c=stats.a_err[:,cc], 
                            vmin=-.2, vmax=.2,
                      edgecolors='k', linewidths=.5, zorder=100, cmap='seismic')
            f.colorbar(sca)
            a.set_title('Amplitude Error (m) | {1} - Obs | {0}'.format(constit[cc], run_name), fontsize=9)
            fn = "ssh_hourly_amp_err_{0}_{1}{2}".format(constit[cc],run_name,file_type)
            f.savefig(os.path.join(dn_out, fn))
            plt.close('all')
            
            # PHASE MAP
            f,a = pu.create_geo_axes(lonbounds, latbounds)
            sca = a.scatter(stats.longitude, stats.latitude, c=stats.g_err[:,cc], 
                            vmin=-15, vmax=15,
                      edgecolors='k', linewidths=.5, zorder=100, cmap='seismic')
            f.colorbar(sca)
            a.set_title('Phase Error (deg) | {1} - Obs | {0}'.format(constit[cc], run_name), fontsize=9)
            fn = "ssh_hourly_pha_err_{0}_{1}{2}".format(constit[cc],run_name,file_type)
            f.savefig(os.path.join(dn_out, fn))
            plt.close('all')
            
            # AMPLITUDE SCATTER
            f,a = pu.scatter_with_fit(stats.a_mod[:,cc], stats.a_obs[:,cc])
            a.set_title('Amplitude Comparison | {0} | {1}'.format(constit[cc], run_name), fontsize=9)
            a.set_xlabel("Model Amplitude (m)")
            a.set_ylabel("Observed Amplitude (m)")
            fn = "ssh_hourly_scatter_amp_{0}_{1}{2}".format(constit[cc], run_name, file_type)
            f.savefig(os.path.join(dn_out, fn))
            plt.close('all')
            
            # PHASE SCATTER
            f,a = pu.scatter_with_fit(stats.g_mod[:,cc], stats.g_obs[:,cc])
            a.set_title('Phase Comparison | {0} | {1}'.format(constit[cc], run_name), fontsize=9)
            a.set_xlabel("Model Phase(deg)")
            a.set_ylabel("Observed Phase (deg)")
            fn = "ssh_hourly_scatter_pha_{0}_{1}{2}".format(constit[cc], run_name, file_type)
            f.savefig(os.path.join(dn_out, fn))
            plt.close('all')
        
        ### CLIM VAR ERROR
        # Jan
        m_ind = [0,3,6,9]
        m_str = ['Jan','Mar','Jul','Oct']
        n_m = len(m_ind)
        for mm in range(0,n_m):
            f,a = pu.create_geo_axes(lonbounds, latbounds)
            sca = a.scatter(stats.longitude, stats.latitude, c=stats.ntr_mod_clim_var[:,m_ind[mm]] - stats.ntr_obs_clim_var[:,m_ind[mm]], 
                            vmin=-.05, vmax=.05,
                      edgecolors='k', linewidths=.5, zorder=100, cmap='seismic')
            f.colorbar(sca)
            a.set_title('Climatological Variance Error | {1} | var({0}) - var(Obs)'.format(run_name, m_str[mm]), fontsize=9)
            fn = "ssh_hourly_clim_var_error_{0}_{1}{2}".format(m_str[mm],run_name,file_type)
            f.savefig(os.path.join(dn_out, fn))
            plt.close('all')
        
        ### FITTED SCATTER PLOTS
        # Plot st.dev comparisons
        f,a = pu.scatter_with_fit(stats.skew_mod, stats.skew_obs)
        a.set_title('Skew Surge Comparison | {0}'.format(run_name), fontsize=9)
        a.set_xlabel("{0} SSH st. dev (m)".format(run_name))
        a.set_ylabel("PSMSL SSH st. dev (m)")
        fn = "ssh_hourly_skew_scatter_{0}{1}".format(run_name, file_type)
        f.savefig(os.path.join(dn_out, fn))
        plt.close('all')
        return
