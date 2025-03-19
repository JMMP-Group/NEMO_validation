"""
module add jaspy
export CONDA_ENV="/home/users/jelt/miniconda3/envs/coast_dev"
source activate $CONDA_ENV
"""

import coast
from coast import general_utils as gu
from coast import plot_util as pu
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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

def analyse_ssh(fn_ext, fn_out, thresholds = np.arange(-.4, 2, 0.1),
                constit_to_save = ['M2','S2','K2','N2','K1','O1','P1','Q1'], 
                semidiurnal_constit = ['M2','S2','K2','N2'],
                diurnal_constit = ['K1','O1','P1','Q1'],
                apply_ntr_filter = True,
		port_id=None ):
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
     port_id              : if None analyse all ports, else (int) analyse selected
    '''
    min_datapoints=4380  # 6month of hrly data

    # Create the object and then inset the netcdf dataset
    if(1): #try:
      ds = xr.open_mfdataset(fn_ext, concat_dim="t_dim",
                                        combine='nested',
                                        #concat_dim="t_dim",
                                        #parallel=True
					)
      ds_ssh = coast.Tidegauge(dataset=ds)
    if(0): #except:
      ds_ssh = coast.Tidegauge(dataset=xr.open_mfdataset(fn_ext, 
					combine='nested', 
					concat_dim="time",
					parallel=True)
				)
    try: ds_ssh.dataset = ds_ssh.dataset.swap_dims({"port":"id_dim"})
    except: pass
    try: ds_ssh.dataset = ds_ssh.dataset.swap_dims({"time":"t_dim"})
    except: pass
    try: ds_ssh.dataset = ds_ssh.dataset.drop_vars("bad_flag")
    except: pass
    
    # Select the target port if specified
    if port_id is not None:
      ds_ssh = ds_ssh.isel(id_dim=[port_id])  # port_id passed as a list so id_dim is not dropped

    # Drop ports that have too few points. Bugged for me due to dask not liking boolean indexing.
    #keep_indices = np.isfinite(ds_ssh.dataset.ssh_obs).sum(dim="t_dim") >= min_datapoints
    #print(keep_indices.values)
    #ds_ssh = ds_ssh.isel(id_dim=keep_indices)
    #print(f"Start with {ds_ssh.dataset.dims['id_dim']} ports")
    #print(f"Keep {keep_indices.values.sum()} ports following check on limited obs")
    #keep_indices = np.isfinite(ds_ssh.dataset.ssh_mod).sum(dim="t_dim") >= min_datapoints
    #ds_ssh = ds_ssh.isel(id_dim=keep_indices)
    #print(f"Keep {keep_indices.values.sum()} ports following check on limited model output")


    # Define Dimension Sizes
    n_port = ds_ssh.dataset.dims['id_dim']
    n_time = ds_ssh.dataset.dims['t_dim']
    n_constit = len(constit_to_save)
    n_thresholds = len(thresholds)
    seasons = ['DJF','JJA','MAM','SON','All']
    n_seasons = len(seasons)
    freq_bands = ['diurnal', 'semidiurnal', 'all']
    n_freq_bands = len(freq_bands)

    # Identify seasons
    month_season_dict = {1:0, 2:0, 3:2, 4:2, 5:2, 6:1,
                         7:1, 8:1, 9:3, 10:3, 11:3, 12:0}
    time = ds_ssh.dataset.time.values
    pd_time = pd.to_datetime(time)
    pd_month = pd_time.month
    pd_season = [month_season_dict[ii] for ii in pd_month]
    
    # NTR dataset
    ds_ntr = xr.Dataset(coords = dict(
                            time = ('t_dim', ds_ssh.dataset.time.values),
                            longitude = ('id_dim', ds_ssh.dataset.longitude.values),
                            latitude = ('id_dim', ds_ssh.dataset.latitude.values)),
                        data_vars = dict(
                            ntr_mod = (['id_dim','t_dim'], np.zeros((n_port, n_time))*np.nan),
                            ntr_obs = (['id_dim','t_dim'], np.zeros((n_port, n_time))*np.nan),
                            ntr_err = (['id_dim','t_dim'], np.zeros((n_port, n_time))*np.nan),
                            ntr_square_err = (['id_dim','t_dim'], np.zeros((n_port, n_time))*np.nan),
                            ntr_abs_err = (['id_dim','t_dim'], np.zeros((n_port, n_time))*np.nan)))
    
    ds_tide = xr.Dataset(coords = dict(
                            time = ('t_dim', ds_ssh.dataset.time.values),
                            longitude = ('id_dim', ds_ssh.dataset.longitude.values),
                            latitude = ('id_dim', ds_ssh.dataset.latitude.values),
                            freq_band = ('freq_band', freq_bands)),
                        data_vars = dict(
                            tide_mod = (['id_dim','freq_band','t_dim'], np.zeros((n_port, n_freq_bands, n_time))*np.nan),
                            tide_obs = (['id_dim','freq_band','t_dim'], np.zeros((n_port, n_freq_bands, n_time))*np.nan),
                            tide_err = (['id_dim','freq_band','t_dim'], np.zeros((n_port, n_freq_bands, n_time))*np.nan),
                            tide_square_err = (['id_dim','freq_band','t_dim'], np.zeros((n_port, n_freq_bands, n_time))*np.nan),
                            tide_abs_err = (['id_dim','freq_band','t_dim'], np.zeros((n_port, n_freq_bands,n_time))*np.nan)))
    
    # ANALYSIS dataset
    ds_stats = xr.Dataset(coords = dict(
                    longitude = ('id_dim', ds_ssh.dataset.longitude.values),
                    latitude = ('id_dim', ds_ssh.dataset.latitude.values),
                    time = ('t_dim', ds_ssh.dataset.time.values),
                    season = ('season', seasons),
                    constituent = ('constituent', constit_to_save),
                    threshold = ('threshold', thresholds)),
               data_vars = dict(
                    a_mod = (['id_dim','constituent'], np.zeros((n_port, n_constit))*np.nan),
                    a_obs = (['id_dim','constituent'], np.zeros((n_port, n_constit))*np.nan),
                    g_mod = (['id_dim','constituent'], np.zeros((n_port, n_constit))*np.nan),
                    g_obs = (['id_dim','constituent'], np.zeros((n_port, n_constit))*np.nan),
                    a_err = (['id_dim','constituent'], np.zeros((n_port, n_constit))*np.nan),
                    g_err = (['id_dim','constituent'], np.zeros((n_port, n_constit))*np.nan),
                    ssh_std_obs = (['id_dim','season'], np.zeros((n_port, n_seasons))*np.nan),
                    ssh_std_mod = (['id_dim','season'], np.zeros((n_port, n_seasons))*np.nan),
                    ssh_std_err = (['id_dim','season'], np.zeros((n_port, n_seasons))*np.nan),
                    ntr_corr = (['id_dim','season'],  np.zeros((n_port, n_seasons))*np.nan),
                    ntr_mae  = (['id_dim','season'], np.zeros((n_port, n_seasons))*np.nan),
                    ntr_me  = (['id_dim','season'], np.zeros((n_port, n_seasons))*np.nan),
                    ntr_rmse = (['id_dim','season'],  np.zeros((n_port, n_seasons))*np.nan),
                    ntr_err_std = (['id_dim','season'],  np.zeros((n_port, n_seasons))*np.nan),
                    ntr_std_obs = (['id_dim','season'],   np.zeros((n_port, n_seasons))*np.nan),
                    ntr_std_mod = (['id_dim','season'],   np.zeros((n_port, n_seasons))*np.nan),
                    thresh_peak_mod  = (['id_dim', 'threshold'], np.zeros((n_port, n_thresholds))),
                    thresh_peak_obs  = (['id_dim', 'threshold'], np.zeros((n_port, n_thresholds))),
                    thresh_time_mod = (['id_dim', 'threshold'], np.zeros((n_port, n_thresholds))),
                    thresh_time_obs = (['id_dim', 'threshold'], np.zeros((n_port, n_thresholds))),
                    thresh_dailymax_mod = (['id_dim', 'threshold'], np.zeros((n_port, n_thresholds))),
                    thresh_dailymax_obs = (['id_dim', 'threshold'], np.zeros((n_port, n_thresholds))),
                    thresh_monthlymax_mod = (['id_dim', 'threshold'], np.zeros((n_port, n_thresholds))),
                    thresh_monthlymax_obs = (['id_dim', 'threshold'], np.zeros((n_port, n_thresholds))))) 
    


    # Commence analysis
    ###################
    tganalysis = coast.TidegaugeAnalysis()

    # Subtract means from all time series
    print(ds_ssh.dataset)
    print(ds_ssh.dataset.mean(dim="t_dim"))
    ds = tganalysis.demean_timeseries(ds_ssh.dataset)

    # Harmonic analysis
    ha_mod = tganalysis.harmonic_analysis_utide(ds.dataset.ssh_mod, min_datapoints=min_datapoints)
    ha_obs = tganalysis.harmonic_analysis_utide(ds.dataset.ssh_obs, min_datapoints=min_datapoints)

    for pp in range(n_port):
        # save constituents (name + amp/pha)
        a_dict_obs = dict(zip(ha_obs[pp].name, ha_obs[pp].A))
        a_dict_mod = dict(zip(ha_mod[pp].name, ha_mod[pp].A))
        g_dict_obs = dict(zip(ha_obs[pp].name, ha_obs[pp].g))
        g_dict_mod = dict(zip(ha_mod[pp].name, ha_mod[pp].g))

        for cc in range(0, len(constit_to_save)):
            if constit_to_save[cc] in ha_obs[pp].name:
                ds_stats['a_mod'][pp, cc] = a_dict_mod[constit_to_save[cc]]
                ds_stats['a_obs'][pp, cc] = a_dict_obs[constit_to_save[cc]]
                ds_stats['a_err'][pp, cc] = a_dict_mod[constit_to_save[cc]] - a_dict_obs[constit_to_save[cc]]
                ds_stats['g_mod'][pp, cc] = g_dict_mod[constit_to_save[cc]]
                ds_stats['g_obs'][pp, cc] = g_dict_obs[constit_to_save[cc]]
                ds_stats['g_err'][pp, cc] = g_dict_mod[constit_to_save[cc]] - g_dict_obs[constit_to_save[cc]]


    # Reconstruct full tidal signal
    tide_mod = tganalysis.reconstruct_tide_utide(ds.dataset.time, ha_mod)
    tide_obs = tganalysis.reconstruct_tide_utide(ds.dataset.time, ha_obs)

    # Reconstruct partial semidiurnal tidal signal
    tide_2_obs = tganalysis.reconstruct_tide_utide(ds.dataset.time, ha_obs,
                                         constit=semidiurnal_constit)
    tide_2_mod = tganalysis.reconstruct_tide_utide(ds.dataset.time, ha_mod,
                                         constit=semidiurnal_constit)

    # Reconstruct partial diurnal tidal signal
    tide_1_obs = tganalysis.reconstruct_tide_utide(ds.dataset.time, ha_obs,
                                         constit=diurnal_constit)
    tide_1_mod = tganalysis.reconstruct_tide_utide(ds.dataset.time, ha_mod,
                                         constit=diurnal_constit)
    
    # Compute differences
    tide_err = tganalysis.difference(tide_mod.dataset, tide_obs.dataset)

    # Store tidal reconstruction data
    for pp in range(n_port):
        # Store full tidal reconstruction
        ds_tide['tide_obs'][pp, -1, :] = tide_obs.isel(id_dim=pp).dataset.reconstructed
        ds_tide['tide_mod'][pp, -1, :] = tide_mod.isel(id_dim=pp).dataset.reconstructed

        ds_tide['tide_err'][pp, -1, :] = tide_err.isel(id_dim=pp).dataset.diff_reconstructed
        ds_tide['tide_abs_err'][pp, -1, :] = tide_err.isel(id_dim=pp).dataset.abs_diff_reconstructed
        ds_tide['tide_square_err'][pp, -1, :] = tide_err.isel(id_dim=pp).dataset.square_diff_reconstructed

        # Store semidiurnal tidal reconstruction
        ds_tide['tide_obs'][pp, 1, :] = tide_2_obs.isel(id_dim=pp).dataset.reconstructed
        ds_tide['tide_mod'][pp, 1, :] = tide_2_mod.isel(id_dim=pp).dataset.reconstructed

        # Store diurnal tidal reconstruction
        ds_tide['tide_obs'][pp, 0, :] = tide_1_obs.isel(id_dim=pp).dataset.reconstructed
        ds_tide['tide_mod'][pp, 0, :] = tide_1_mod.isel(id_dim=pp).dataset.reconstructed

    # NTR (non tidal residual)
    ntr_mod = tganalysis.calculate_non_tidal_residuals(ds.dataset.ssh_mod, tide_mod.dataset.reconstructed,
                                                       apply_filter=False)
    ntr_obs = tganalysis.calculate_non_tidal_residuals(ds.dataset.ssh_mod, tide_obs.dataset.reconstructed,
                                                       apply_filter=True, window_length=10, polyorder=2)
    # Calculate errors
    ntr_diff = tganalysis.difference(ntr_mod.dataset, ntr_obs.dataset)
    ds_ntr['ntr_obs'] = ntr_obs.dataset.ntr
    ds_ntr['ntr_mod'] = ntr_mod.dataset.ntr
    ds_ntr['ntr_err'] = ntr_diff.dataset.diff_ntr
    ds_ntr['ntr_abs_err'] = ntr_diff.dataset.abs_diff_ntr
    ds_ntr['ntr_square_err'] = ntr_diff.dataset.square_diff_ntr

    # Make masked arrays for seasonal correlation calculation
    ntr_corr = xr.corr(ntr_obs.dataset.ntr, ntr_mod.dataset.ntr, dim="t_dim")
    ds_stats['ntr_corr'][:,-1] = ntr_corr

    #  TWL (total water level) standard deviations
    ssh_std = ds_ssh.dataset.std(dim='t_dim', skipna=True)
    ds_stats['ssh_std_mod'][:, 4] = ssh_std.ssh_mod
    ds_stats['ssh_std_obs'][:, 4] = ssh_std.ssh_obs
    ds_stats['ssh_std_err'][:, 4] = ssh_std.ssh_mod - ssh_std.ssh_obs

    # Threshold statistics, on NTR (Slow)
    thresh_mod = tganalysis.threshold_statistics(ntr_mod.dataset, thresholds=thresholds)
    thresh_obs = tganalysis.threshold_statistics(ntr_obs.dataset, thresholds=thresholds)

    ds_stats['thresh_peak_mod'] = thresh_mod.peak_count_ntr
    ds_stats['thresh_peak_obs'] = thresh_obs.peak_count_ntr
    ds_stats['thresh_time_mod'] = thresh_mod.time_over_threshold_ntr
    ds_stats['thresh_time_obs'] = thresh_obs.time_over_threshold_ntr
    ds_stats['thresh_dailymax_mod'] = thresh_mod.dailymax_count_ntr
    ds_stats['thresh_dailymax_obs'] = thresh_obs.dailymax_count_ntr
    ds_stats['thresh_monthlymax_mod'] = thresh_mod.monthlymax_count_ntr
    ds_stats['thresh_monthlymax_obs'] = thresh_obs.monthlymax_count_ntr

    # Merge datasets and write to file
    ds_merge = xr.merge( (ds_stats, ds_ntr, ds_tide))
    write_ds_to_file(ds_merge, fn_out)
    
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
    min_datapoints=4380 # 6 months of hrly  
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
    print(f"fn nemo:{fn_nemo_data}")    
    nemo = coast.Gridded(fn_nemo_data, fn_nemo_domain, config=fn_nemo_cfg,
                              multiple=True)

    # Rename depth_0 to be depth
    nemo.dataset = nemo.dataset.rename({"depth_0": "depth"})
    # Create a landmask array and put it into the nemo object.
    # Here, using the bottom_level == 0 variable from the domain file is enough.
    nemo.dataset["landmask"] = nemo.dataset.bottom_level == 0

    # Read OBS data from tg file
    obs = coast.Tidegauge(dataset=xr.open_dataset(fn_obs))
    obs.dataset = obs.dataset.set_coords("time")
    obs.dataset = obs.dataset.swap_dims({"port": "id_dim"})
    obs.dataset = obs.dataset.swap_dims({"time": "t_dim"})  # expected dim in COAsT, though some xarray func requires time

    # Cut the OBS data down the to the same time window as the model
    start_date = nemo.dataset.time.isel(t_dim=0).values
    end_date = nemo.dataset.time.isel(t_dim=-1).values
    obs = obs.time_slice( date0=start_date, date1=end_date )

    # Then do the interpolation
    model_timeseries = obs.obs_operator(nemo, time_interp='linear')

    # Drop points that are too far away from observation
    keep_indices = model_timeseries.dataset.interp_dist <= dist_omit
    model_timeseries = model_timeseries.isel(id_dim=keep_indices)
    # Also drop the unwanted observational timeseries
    obs = obs.isel(id_dim=keep_indices)

    tganalysis = coast.TidegaugeAnalysis()
    # This routine searches for missing values in each dataset and applies them
    # equally to each corresponding dataset
    obs_new, model_new = tganalysis.match_missing_values(obs.dataset.ssh, model_timeseries.dataset.ssh)

    # Take shared coords from obs. Merge ssh values and save
    if np.abs(model_new.dataset.time.values - obs_new.dataset.time.values).max() == np.timedelta64(0, 'ns'):
        #model_new.dataset['longitude'] = obs_new.dataset.longitude
        #model_new.dataset['latitude'] = obs_new.dataset.latitude
        #model_new.dataset['time'] = obs_new.dataset.time
        #ds_extract = xr.merge((obs_new.dataset, model_new.dataset))
        ds_extract = obs_new
        ds_extract.dataset = ds_extract.dataset.rename_vars({"ssh": "ssh_obs"})
        ds_extract.dataset['ssh_mod'] = model_new.dataset.ssh
        write_ds_to_file(ds_extract.dataset, fn_out)
        print(f"Extracted dataset saved to: {fn_out}")
    else:
        print(f"Need to fix issue with common time axis. Expected identical time values")
        print(f"Diff: {np.abs(model_new.dataset.time.values - obs_new.dataset.time.values).max()}")
        print(f"Extracted dataset not saved")


    if(0):
        print()


        # Subset obs by NEMO domain
        lonmax = np.nanmax(nemo.longitude)
        lonmin = np.nanmin(nemo.longitude)
        latmax = np.nanmax(nemo.latitude)
        latmin = np.nanmin(nemo.latitude)
        ind = gu.subset_indices_lonlat_box(obs.longitude, obs.latitude,
                                           lonmin, lonmax, latmin, latmax)
        obs = obs.isel(port=ind[0])  # only the obs in the box

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

class TaylorTide():
    def __init__(self,
                r_obs = None,
                rms_amp_max = 0.61,
                rms_amp_contours = [0.2, 0.4, 0.6],
                rms_err_contours=[0.2, 0.4, 0.6],
                cos_theta_lines = [0.3, 0.6, 0.9],
                theta_lines = [0, 15, 30, 45, 60],
                theta_lines_flag = True, # theta or cos(theta) construction lines
                err_contour_flag = True,
                ):

        self.fig, self.ax = self.plot_frame(r_obs, rms_amp_max, rms_amp_contours,
                                            rms_err_contours, cos_theta_lines, theta_lines,
                                            theta_lines_flag, err_contour_flag)

    # This custom formatter removes trailing zeros, e.g. "1.0" becomes "1", and
    def rms_fmt(self,x):
        s = f"{x:.2f}"
        if s.endswith("00"):  # ends with two zeros
            s = f"{x:.0f}"
        elif s.endswith("0"): # end with ONLY one zero
            s = f"{x:.1f}"
        return rf"{s}"

    def plot_frame(self, r_obs, rms_amp_max, rms_amp_contours, rms_err_contours,
                   cos_theta_lines, theta_lines,
                   theta_lines_flag, err_contour_flag):

        #fig = plt.figure(figsize=(2, 2)) # to make a thumbnail schematic
        fig = plt.figure()

        theta = np.arange(0, np.pi + np.pi/100, np.pi/100)

        # define meshgrid for contour plots
        x = np.arange(0,rms_amp_max, rms_amp_max/100)
        y = np.arange(0,rms_amp_max, rms_amp_max/100)
        X,Y = np.meshgrid(x,y)

        # setting the axis limits in [left, bottom, width, height]
        rect = [0.1, 0.1, 0.8, 0.8]

        # the cartesian axis:
        ax = fig.add_axes(rect, frameon=False)
        #ax =fig.add_subplot(111)

        # RMS amplitude arc contours
        if rms_amp_contours != []:
            Camp = ax.contour( X,Y,np.sqrt(X**2 + Y**2), levels=rms_amp_contours, colors='grey', linestyles='dotted')
            ax.clabel(Camp, Camp.levels, inline=True, fmt=self.rms_fmt, fontsize=10)

        # Obs point, arc through obs, RMS error arcs
        if r_obs is not None:
            ax.scatter(r_obs, 0, s=30, color='blue', clip_on=False)  # obs point

            if err_contour_flag is True:
                ax.plot(r_obs * np.cos(theta), r_obs * np.sin(theta), '-', color='blue')  # arc through obs

                # RMS error arc from obs as origin
                mask = X**2 + Y**2 > rms_amp_max**2
                C = np.ma.masked_where(mask, np.sqrt((X-r_obs)**2 + Y**2))
                Cerr = ax.contour(X, Y, C, levels=rms_err_contours, colors='grey',
                                  linestyles='dashed')
                ax.clabel(Cerr, Cerr.levels, inline=True, fmt=self.rms_fmt, fontsize=10)

                # Add text contour label. THIS WILL BE A PROBLEM IF MORE THAN 3 LEVELS ARE DEFINED
                fmt = {}
                strs = ['', '', 'RMS error']
                for l, s in zip(rms_err_contours, strs):
                    fmt[l] = s
                ax.clabel(Cerr, Cerr.levels, inline=True, fmt=fmt, fontsize=10)


        # Bounding lines - black
        ax.plot( rms_amp_max*np.cos(theta), rms_amp_max*np.sin(theta), '-', color='black')
        # Bounding x=0 and y=0 lines - black
        ax.plot( [0,rms_amp_max], [0,0], '-', color='black')
        ax.plot( [0,0], [0,rms_amp_max], '-', color='black')


        # Cos theta / correlation lines
        r = np.arange(0, rms_amp_max, rms_amp_max/50)
        if theta_lines_flag == False:
            for cos_theta in cos_theta_lines:
                ax.plot( r*cos_theta, r*np.sqrt(1 - cos_theta**2), ':', color='grey')
                ax.text( rms_amp_max*cos_theta, rms_amp_max*np.sqrt(1 - cos_theta**2), str(cos_theta), color='k' )
            ax.text( rms_amp_max*1/np.sqrt(2), rms_amp_max*1/np.sqrt(2), "Correlation", rotation=-45, color='k')
        else:
            for ang in theta_lines:
                ang_rad = ang*np.pi/180.
                ax.plot( r*np.cos(ang_rad), r*np.sin(ang_rad), ':', color='grey')
                ax.text( rms_amp_max*np.cos(ang_rad), rms_amp_max*np.sin(ang_rad), str(ang)+"$^o$", color='k' )
            ax.text( rms_amp_max*1/np.sqrt(2), rms_amp_max*1/np.sqrt(2), "phase error", rotation=-45, color='k')


        # axis limits and axis labels
        ax.set_xlim([0,rms_amp_max])
        ax.set_ylim([0,rms_amp_max])
        ax.set_xlabel('RMS amplitude (m)')
        ax.set_ylabel('RMS amplitude (m)')

        ax.set_aspect(1)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)

        if(0):
            # the polar axis:
            ax_polar = fig.add_axes([0.1, 0.0, 0.8, 0.8], polar=True, frameon=False)
            #ax_polar = fig.add_axes(ax.get_position(original=True), polar=True, frameon=False)
            ## [left, bottom, width, height]
            #ax_polar.set_position([0.1, 0.0, 0.8, 0.8])

            # the polar plot
            ax_polar.plot(r, r, color='r', linewidth=3)
            ax_polar.set_rmax(2.0)
            ax_polar.grid(True)
            ax_polar.set_thetamin(0)
            ax_polar.set_thetamax(54)
            ax_polar.set_xlabel('RMS amplitude (m)')
            ax_polar.set_ylabel('RMS amplitude (m)')


        return fig, ax

class plot_taylor_tide():

    def __init__(self, fn_ssh_hourly_stats, dn_out, run_name, file_type='.png'):
# define colors
        self.c_co7 = 'lightblue' #[46/256., 153/256., 195/256.]
        self.c_co9 = 'green'
        #c_FES = 'y'
        self.c_obs = 'blue'
        self.plot_fn(       fn_ssh_hourly_stats, dn_out, run_name, file_type='.png')

    def pearson_correl_coef(self, z1obs, z2obs, z1mod, z2mod):
       """
       sum( <zobs_i, zmod_i> )/ ( rms(zobs) rms(zmod) )
        for vectors zobs_i = (z1obs_i, z2obs_i))
        and zero mean zobs and zmod
       """
       inner_product = np.nanmean( z1obs*z1mod + z2obs*z2mod )
       return inner_product / (self.rms_abs_amp(z1obs,z2obs) * self.rms_abs_amp(z1mod,z2mod))

    def rms_abs_amp(self, z1, z2):
      """
      root mean square vector amplitude
      z1=a.cos(theta), z2=a.sin(theta)
      """
      return np.sqrt(np.nanmean( z1**2 + z2**2 ))

    def rms_abs_error(self, z1obs, z2obs, z1mod, z2mod):
      """
      root mean square absolute error = sqrt( sum |zmod - zobs|^2 )
      """
      return self.rms_abs_amp( z1obs - z1mod, z2obs - z2mod )

    def plot_fn(self, fn_ssh_hourly_stats, dn_out, run_name, file_type='.png'):
        
        if type(fn_ssh_hourly_stats) is not list:
            fn_ssh_hourly_stats = [fn_ssh_hourly_stats]
            run_name = [run_name]

        if len(run_name) != len(fn_ssh_hourly_stats):
            print('Error: run_name length must be equal to fn_ssh_hourly_stats length.')
            return

        n_cfgs = len(fn_ssh_hourly_stats)

        stats = [xr.open_dataset(fn ,engine="netcdf4") for fn in fn_ssh_hourly_stats]

        #stats = xr.open_dataset(fn_ssh_hourly_stats, engine="netcdf4")

        ## load data
        #a_obs = stats[0].a_obs  # [nobs, constit]
        #g_obs = stats[0].g_obs  # [nobs, constit]
        #a_mod = stats[0].a_mod  # [nobs, constit]
        #g_mod = stats[0].g_mod  # [nobs, constit]
       
        # Loop over harmonic species.
        # Plot Taylor Tide diag of model and obs for each harmonic. Overlay on two (deep/shallow) plots as trees

        nsim=3 # number of simulations + obs.  labels = ["co9", "co7", "obs"]

        #for constit_family_list in [["M2", "S2", "N2", "K2"]]:
        for constit_family_list in [["M2", "S2", "N2", "K2"], ["O1", "Q1", "P1"]]:
          if "M2" in constit_family_list: family_str = "semi-diurnal"
          if "O1" in constit_family_list: family_str = "diurnal"
          R = np.zeros((nsim, len(constit_family_list)))
          rms_amp = np.zeros((nsim, len(constit_family_list)))
          rms_err = np.zeros((nsim, len(constit_family_list)))
          label = {}
          for count, constit in enumerate(constit_family_list):
            try:
                del II, z1obs, z2obs, z1mod, z2mod
            except:
                pass

            cc = np.argwhere( (stats[0].constituent == constit).values)[0][0] # find index of constit
            z1obs = stats[0].a_obs[:,cc]*np.cos(np.deg2rad(stats[0].g_obs[:,cc]))
            z2obs = stats[0].a_obs[:,cc]*np.sin(np.deg2rad(stats[0].g_obs[:,cc])) 

            R[0,count] = 1
            rms_amp[0,count] = self.rms_abs_amp(z1obs, z2obs)
            rms_err[0,count] = 0
            label[0,count] = 'obs:'+constit


            # CO7_AMM15
            try:
                del z1mod, z2mod
            except:
                pass
            cc = np.argwhere( (stats[0].constituent == constit).values)[0][0] # find index of constit
            z1mod = stats[0].a_mod[:,cc]*np.cos(np.deg2rad(stats[0].g_mod[:,cc]))
            z2mod = stats[0].a_mod[:,cc]*np.sin(np.deg2rad(stats[0].g_mod[:,cc]))

            R[1,count] = self.pearson_correl_coef(z1obs, z2obs, z1mod, z2mod)
            rms_amp[1,count] = self.rms_abs_amp(z1mod, z2mod)
            rms_err[1,count] = self.rms_abs_error(z1obs, z2obs, z1mod, z2mod)
            label[1,count] = 'CO7_AMM15:'+constit

            if n_cfgs == 2:
              # CO9_AMM15
              try:
                del z1mod, z2mod
              except:
                pass
              cc = np.argwhere( (stats[1].constituent == constit).values)[0][0] # find index of constit
              z1mod = stats[1].a_mod[:,cc]*np.cos(np.deg2rad(stats[1].g_mod[:,cc]))
              z2mod = stats[1].a_mod[:,cc]*np.sin(np.deg2rad(stats[1].g_mod[:,cc]))

              R[2,count] = self.pearson_correl_coef(z1obs, z2obs, z1mod, z2mod)
              rms_amp[2,count] = self.rms_abs_amp(z1mod, z2mod)
              rms_err[2,count] = self.rms_abs_error(z1obs, z2obs, z1mod, z2mod)
              label[2,count] = 'CO9_AMM15:'+constit


          print(label)

          ## Check cosine rule consistency. Output data.
          for count, constit in enumerate(constit_family_list):
          #for i in range(len(constit_family_list)):
                B = rms_amp[0,count]
                A = rms_amp[1:nsim,count]
                C = rms_err[1:nsim,count]
                costheta = R[1:nsim,count]

                #print(constit)
                #print(B)
                #print(A)
                #print(C)
                #print(costheta)
                print("Check cosine rule consistency")
                for j in range(0,nsim-1): # model runs: fes, gs1p1, gs1p2, tdiss
                    print(
                        f"{label[1+j,count]}: sqrt(A^2+B^2-2ABcos(theta))={np.sqrt(A[j] ** 2 + B ** 2 - 2 * A[j] * B * costheta[j])}. C={C[j]}")

                print("Output table of data")
                for j in range(0, nsim - 1):  # model runs: fes, gs1p1, gs1p2, tdiss
                    print(
                        f"{label[1 + j, count]:11s}: (A, B, C, theta)={A[j]:.3f}, {B:.3f}, {C[j]:.3f}, {np.arccos(costheta[j]) * 180 / np.pi:.1f}")
                del B, A, C, costheta



          # Create TaylorTide plot template
          if family_str == "semi-diurnal":
                fig_amp_max = 2.5
                fig_err_contours = [0.2, 0.4, 0.6]

          elif family_str == "diurnal":
                fig_amp_max = 0.11
                fig_err_contours = [0.01, 0.02, 0.03]



          tt = TaylorTide(
                r_obs=rms_amp[0,0],
                rms_amp_max=fig_amp_max,
                rms_amp_contours=[],
                #rms_amp_contours=[0.2, 0.4, 0.6],
                rms_err_contours=fig_err_contours,
                cos_theta_lines=[0.3, 0.6, 0.9],
                err_contour_flag=True,
            )

          ## Loop over constituents to add data
          for i in range(len(constit_family_list)):
                # Add data to axes
                tt.ax.scatter(rms_amp[1:nsim,i] * R[1:nsim,i],
                          rms_amp[1:nsim,i] * np.sqrt(1 - R[1:nsim,i] ** 2),
                          s=60, c=[ self.c_co7, self.c_co9 ], zorder=10, clip_on=False)
                            #s = 20, c = ['b', 'g', 'r', 'y'], zorder = 10, clip_on = False)
                # Add vectors between points and obs
                tt.ax.plot([np.repeat(rms_amp[0,i],nsim-1), rms_amp[1:nsim,i] * R[1:nsim,i]],
                       [np.zeros(nsim-1), rms_amp[1:nsim,i] * np.sqrt(1 - R[1:nsim,i] ** 2)])
                if nsim != 3: print('Colours and lines not as expected here')
                #tt.ax.lines[-4].set_color('k')
                #tt.ax.lines[-3].set_color(c_ZPS_TIDE) #'g')
                tt.ax.lines[-2].set_color(self.c_co7) #'r')
                tt.ax.lines[-1].set_color(self.c_co9) #'y')

                # Place constituent labels
                #tt.ax.text( rms_amp[0,i], -0.036*fig_amp_max, constit_family_list[i], rotation=0, color='b')
                xpos = rms_amp[0,i]
                ypos = 0.036*fig_amp_max
                if constit_family_list[i]=="M2": xpos = xpos - 0.20
                if constit_family_list[i]=="S2": xpos = xpos + 0.00
                if constit_family_list[i]=="K2": xpos = xpos - 0.03; ypos = ypos + 0.02
                if constit_family_list[i]=="N2": xpos = xpos - 0.02; ypos = ypos + 0.04
                if constit_family_list[i]=="Q1": xpos = xpos - 0.012
                if constit_family_list[i]=="O1": xpos = xpos + 0.002
                if constit_family_list[i]=="P1": xpos = xpos + 0.002
                tt.ax.text( xpos, ypos, constit_family_list[i], fontsize=12,  rotation=0, color='b')

                # Draw obs dot (smaller)
                tt.ax.scatter(rms_amp[0,i], 0, s=40, color=self.c_obs, zorder=20, clip_on=False)  # obs point

          # Add indicative timing labels
          if family_str == "diurnal":
                tt.ax.text( fig_amp_max*np.cos(np.pi/6), fig_amp_max*np.sin(np.pi/6), "      (2hr)", color='k')
                tt.ax.text( fig_amp_max*np.cos(np.pi/24), fig_amp_max*np.sin(np.pi/24), " (30min)", color='k')
          elif family_str == "semi-diurnal":
                tt.ax.text( fig_amp_max*np.cos(np.pi/6), fig_amp_max*np.sin(np.pi/6),  "      (1hr)", color='k')
                tt.ax.text( fig_amp_max*np.cos(np.pi/24), fig_amp_max*np.sin(np.pi/24)," (15min)", color='k')

          # manual legend (once)
          if(1): # family_str == "semi-diurnal":
                #colors = [ 'green', 'red', 'yellow', 'blue']
                colors = [ self.c_co7, self.c_co9, self.c_obs]
                sizes = [ 6, 6, 4]
                lines = [Line2D([], [], color=colors[ci], markersize=sizes[ci], marker='o', linestyle='None') for ci in range(len(colors))]
                labels = [ "CO7_AMM15", "CO9_AMM15", "observations"]
                plt.legend(lines, labels, loc='upper left', framealpha=1)

          #plt.title(subset + ":" + family_str)
          plt.title(family_str ) 

          fn = "Taylor_ssh_hourly_harmonic_tree_"+family_str+file_type
          plt.savefig(os.path.join(dn_out, fn))
          plt.close('all')





class plot_single_cfg():
    
    def __init__(self, fn_ssh_hourly_stats, dn_out, run_name, file_type='.png'):
    
        stats = xr.open_dataset(fn_ssh_hourly_stats, engine="netcdf4")
        
        lonmax = np.max((stats.longitude + 180)%360 - 180)  # centre around lon=0
        lonmin = np.min((stats.longitude + 180)%360 - 180)  # centre around lon=0
        latmax = np.max(stats.latitude)
        latmin = np.min(stats.latitude)
        lonbounds = [lonmin-4, lonmax+4]
        latbounds = [latmin-4, latmax+4]
        
        ### HARMONIC SUMMARY
        f,a = plt.subplots()
        sca = a.pcolormesh( stats.a_obs )
        f.colorbar(sca)
        a.set_title('Harmonic amplitudes for Tide Gauge station observations', fontsize=9)
        a.set_xticks(np.arange(len(stats['constituent']))+0.5)
        a.set_xticklabels(stats['constituent'].values)
        a.set_ylabel('port')
        fn = "ssh_hourly_harmonic_amp_obs"+run_name+file_type
        f.savefig(os.path.join(dn_out, fn))
        plt.close('all')

        ### GEOGRAPHICAL SCATTER PLOTS
        # Plot correlations
        f,a = pu.create_geo_axes(lonbounds, latbounds)
        sca = a.scatter(stats.longitude, stats.latitude, c=stats.ntr_corr.sel(season="All"), 
                        #vmin=.85, vmax=1,
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
        
        # Plot ntr_err_std
        f,a = pu.create_geo_axes(lonbounds, latbounds)
        sca = a.scatter(stats.longitude, stats.latitude, c=stats.ntr_err_std.sel(season="All"), 
                        vmin=-.15, vmax=.15,
                  edgecolors='k', linewidths=.5, zorder=100, cmap='seismic')
        f.colorbar(sca)
        a.set_title('Hourly NTR St. Dev Error (m) | {0}'.format(run_name), fontsize=9)
        fn = "ssh_hourly_ntr_err_std_"+run_name+file_type
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
            for i in range(len(stats.g_mod[:,cc])):
              if stats.g_obs[i,cc] - stats.g_mod[i,cc] > 250:
                stats.g_mod[i,cc] = stats.g_mod[i,cc] + 360
              if stats.g_mod[i,cc] - stats.g_obs[i,cc] > 300:
                stats.g_mod[i,cc] = stats.g_mod[i,cc] - 360
            f,a = pu.scatter_with_fit(stats.g_mod[:,cc], stats.g_obs[:,cc])
            a.set_title('Phase Comparison | {0} | {1}'.format(constit[cc], run_name), fontsize=9)
            a.set_xlabel("Model Phase(deg)")
            a.set_ylabel("Observed Phase (deg)")
            fn = "ssh_hourly_scatter_pha_{0}_{1}{2}".format(constit[cc], run_name, file_type)
            f.savefig(os.path.join(dn_out, fn))
            plt.close('all')
        
        # THRESHOLD PLOTS
        prop = stats.thresh_peak_mod/stats.thresh_peak_obs
        plot_thresh = [5, 10, 15]
        for pp in range(0, len(plot_thresh)):
            prop_tmp = prop.isel(threshold=plot_thresh[pp])
            f,a = pu.create_geo_axes(lonbounds, latbounds)
            sca = a.scatter(stats.longitude, stats.latitude, c=prop_tmp, 
                            #vmin=0, vmax=2,
                      edgecolors='k', linewidths=.5, zorder=100, cmap='seismic')
            f.colorbar(sca)
            a.set_title('Number of model NTR peaks as a proportion of observed NTR peaks \n Threshold = {0}m | {1}'.format(prop_tmp.threshold.values, run_name), fontsize=9)
            fn = "thresh_freq_{0}_{1}{2}".format(prop_tmp.threshold.values, run_name, file_type)
            f.savefig(os.path.join(dn_out, fn))
            plt.close('all')
            
        prop = stats.thresh_time_mod/stats.thresh_time_obs
        plot_thresh = [5, 10, 15]
        for pp in range(0, len(plot_thresh)):
            prop_tmp = prop.isel(threshold=plot_thresh[pp])
            f,a = pu.create_geo_axes(lonbounds, latbounds)
            sca = a.scatter(stats.longitude, stats.latitude, c=prop_tmp, 
                            #vmin=0, vmax=2,
                      edgecolors='k', linewidths=.5, zorder=100, cmap='seismic')
            f.colorbar(sca)
            a.set_title('Model hours over threshold as a proportion of observed time \n Threshold = {0}m | {1}'.format(prop_tmp.threshold.values, run_name), fontsize=9)
            fn = "thresh_int_{0}_{1}{2}".format(prop_tmp.threshold.values, run_name, file_type)
            f.savefig(os.path.join(dn_out, fn))
            plt.close('all')
        
        
        ### THRESHOLD LINE PLOTS
        f,a = plt.subplots() 
        sca = a.plot(stats.threshold, stats.thresh_time_mod.mean(dim="id_dim"), label=run_name)
        sca = a.plot(stats.threshold, stats.thresh_time_obs.mean(dim="id_dim"), label="obs")
        a.set_title('Duration over threshold as a proportion of total duration | {0}'.format(run_name), fontsize=9)
        a.set_xlim([0,1.5])
        plt.xlabel("Threshold (m)")
        plt.legend()
        fn = "thresh_time_{0}{1}".format(run_name, file_type)
        f.savefig(os.path.join(dn_out, fn))
        plt.show()
        plt.close('all')


        f,a = plt.subplots() 
        sca = a.plot(stats.threshold, stats.thresh_peak_mod.mean(dim="id_dim"), label=run_name)
        sca = a.plot(stats.threshold, stats.thresh_peak_obs.mean(dim="id_dim"), label="obs")
        a.set_title('Peak count over threshold as a proportion of total peaks | {0}'.format(run_name), fontsize=9)
        a.set_xlim([0,1.5])
        plt.xlabel("Threshold (m)")
        plt.legend()
        fn = "thresh_int_{0}{1}".format(run_name, file_type)
        f.savefig(os.path.join(dn_out, fn))
        plt.show()
        plt.close('all')


        
        if(0):
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
        else:
          print(f"Not plotted scatter of clim ntr errors")

        if(0):
          ### FITTED SCATTER PLOTS
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
        else:
          print(f"Not plotted scatter skew comparision")        
        
class plot_stats_ssh_hourly_multi_cfg():
    """ Overlay plots from multi cfgs"""
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
    """
    Plot differences between two configs. 
    Convention: 2nd config is the reference
     """
    def __init__(self, fn_ssh_hourly_stats1, fn_ssh_hourly_stats2, dn_out, 
                 run_name1, run_name2, file_type='.png'):
    
        stats1 = xr.open_dataset(fn_ssh_hourly_stats1)
        stats2 = xr.open_dataset(fn_ssh_hourly_stats2)
       
        lonmax = np.max((stats1.longitude + 180)%360 - 180)  # centre around lon=0
        lonmin = np.min((stats1.longitude + 180)%360 - 180)  # centre around lon=0
        latmax = np.max(stats1.latitude)
        latmin = np.min(stats1.latitude)
        lonbounds = [lonmin-4, lonmax+4]
        latbounds = [latmin-4, latmax+4]
 
        ### GEOGRAPHICAL SCATTER PLOTS
        # Plot difference in NTR correlations
        f,a = pu.create_geo_axes(lonbounds, latbounds)
        sca = a.scatter(stats1.longitude, stats1.latitude, c=(stats1.ntr_corr-stats2.ntr_corr).sel(season="All"),
                        vmin=-.1, vmax=.1,
                  edgecolors='k', linewidths=.5, zorder=100, cmap='RdYlGn')
        f.colorbar(sca)
        a.set_title('Hourly NTR Correlations with Tide Gauge difference | corr({0}) - corr({1})'.format(run_name1,run_name2), fontsize=9)
        fn = "ssh_hourly_ntr_correlations_difference"+run_name1+"-"+run_name2+file_type
        f.savefig(os.path.join(dn_out, fn))
        plt.close('all')

        # Plot ssh_std_err = std(mod) - std(obs)
        f,a = pu.create_geo_axes(lonbounds, latbounds)
        sca = a.scatter(stats1.longitude, stats1.latitude, c=(np.abs(stats2.ssh_std_err)-np.abs(stats1.ssh_std_err)).sel(season="All"),
                        vmin=-.15, vmax=.15,
                  edgecolors='k', linewidths=.5, zorder=100, cmap='RdYlGn')
        f.colorbar(sca)
        a.set_title('Improvement of hourly SSH St. Dev Error (m) | abs({1})-abs({0})'.format(run_name1, run_name2), fontsize=9)
        fn = "ssh_hourly_stderr_diffabs"+run_name2+"-"+run_name1+file_type
        f.savefig(os.path.join(dn_out, fn))
        plt.close('all')

        
        # Plot difference in ntr_err_std
        f,a = pu.create_geo_axes(lonbounds, latbounds)
        sca = a.scatter(stats1.longitude, stats1.latitude, c=(stats2.ntr_err_std-stats1.ntr_err_std).sel(season="All"),
                        #vmin=-.15, vmax=.15,
                  edgecolors='k', linewidths=.5, zorder=100, cmap='RdYlGn')
        f.colorbar(sca)
        a.set_title('Hourly NTR St. Dev Error (m) difference | {1} - {0}'.format(run_name1,run_name2), fontsize=9)
        fn = "ssh_hourly_ntr_stderr_difference"+run_name2+"-"+run_name1+file_type
        f.savefig(os.path.join(dn_out, fn))
        plt.close('all')
        
        # Plot difference in mae
        f,a = pu.create_geo_axes(lonbounds, latbounds)
        sca = a.scatter(stats1.longitude, stats1.latitude, c=(stats2.ntr_mae-stats1.ntr_mae).sel(season="All"),
                        #vmin=-.05, vmax=.05,
                  edgecolors='k', linewidths=.5, zorder=100)
        f.colorbar(sca)
        a.set_title('Hourly NTR MAE (m) difference | {1} - {0}'.format(run_name1,run_name2), fontsize=9)
        fn = "ssh_hourly_ntr_mae_difference"+run_name2+"-"+run_name1+file_type
        f.savefig(os.path.join(dn_out, fn))
        plt.close('all')
        
        ### HARMONIC ERRORS 
        constit = stats1.constituent.values
        n_constit = len(constit)
        
        for cc in range(0,n_constit):
            # AMPLITUDE MAP DIFF of ABS; (abs_err2 - abs_err1)
            f,a = pu.create_geo_axes(lonbounds, latbounds)
            sca = a.scatter(stats1.longitude, stats1.latitude, c=np.abs(stats2.a_err[:,cc])-np.abs(stats1.a_err[:,cc]),
                            vmin=-.05, vmax=.05,
                      edgecolors='k', linewidths=.5, zorder=100, cmap='RdYlGn')
            f.colorbar(sca)
            a.set_title('Amplitude Error (m). Abs diff | {2} - {1} | {0}'.format(constit[cc], run_name1, run_name2), fontsize=9)
            fn = "ssh_hourly_amp_err_diffabs_{0}_{2}-{1}{3}".format(constit[cc], run_name1, run_name2, file_type)
            f.savefig(os.path.join(dn_out, fn))
            plt.close('all')
            
            if(0):
                # PHASE MAP
                f,a = pu.create_geo_axes(lonbounds, latbounds)
                sca = a.scatter(stats1.longitude, stats1.latitude, c=stats.g_err[:,cc],
                                vmin=-15, vmax=15,
                          edgecolors='k', linewidths=.5, zorder=100, cmap='seismic')
                f.colorbar(sca)
                a.set_title('Phase Error (deg) | {1} - Obs | {0}'.format(constit[cc], run_name), fontsize=9)
                fn = "ssh_hourly_pha_err_{0}_{1}{2}".format(constit[cc],run_name,file_type)
                f.savefig(os.path.join(dn_out, fn))
                plt.close('all')
            else:
                if cc == 0: print(f"Diff of phase maps not plotted")

            if(0):
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
            else:
                if cc == 0: print(f"Diff of climatological NTR not plotted")
        
        if(0):
            ### FITTED SCATTER PLOTS
            f,a = pu.scatter_with_fit(stats.skew_mod, stats.skew_obs)
            a.set_title('Skew Surge Comparison | {0}'.format(run_name), fontsize=9)
            a.set_xlabel("{0} SSH st. dev (m)".format(run_name))
            a.set_ylabel("PSMSL SSH st. dev (m)")
            fn = "ssh_hourly_skew_scatter_{0}{1}".format(run_name, file_type)
            f.savefig(os.path.join(dn_out, fn))
            plt.close('all')
        else:
            print(f"Diff of scatter skew not plotted")
        return
