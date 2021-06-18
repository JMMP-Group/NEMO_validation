import xarray as xr
import xarray.ufuncs as uf
import numpy as np
import pandas as pd

import sys
# UNCOMMENT THIS LINE IF USING DEVELOPMENT VERSION OF COAST
sys.path.append('/home/users/dbyrne/code/COAsT')
import coast
import coast.general_utils as gu
import coast.crps_util as cu
import os
import os.path
from datetime import datetime, timedelta

def get_season_index(dt):
    ''' return array of season indices for a given list of datetimes '''
    month_season_dict = {1:1, 2:1, 3:2, 4:2, 5:2, 6:3,
    7:3, 8:3, 9:4, 10:4, 11:4, 12:1}
    dt = pd.to_datetime(dt)
    month_index = dt.month
    season_index = [month_season_dict[mm] for mm in month_index]
    return season_index

def write_ds_to_file(ds, fn, **kwargs):
    ''' Quick routine for writing to file and checking if file already exists '''
    if os.path.exists(fn):
        os.remove(fn)
    ds.to_netcdf(fn, **kwargs)
    
    
def analyse_regional(fn_stats, fn_nemo_domain, fn_out, 
                       regional_masks=[], region_names=[],
                       start_date = None, end_date = None):
                       
    '''
    A routine for averaging the hourly analysis into regional subdomains.
    This routine will also average into seasons as well as over all time.
    The resulting output dataset will have dimension (region, season)
    
    INPUTS
     fn_stats (str)        : Absolute path to analysis output from analyse_ts_hourly_en4()
     fn_nemo_domain (str)  : Absolute path to NEMO domain_cfg.nc
     fn_out (str)          : Absolute path to desired output file
     regional_masks (list) : List of 2D boolean arrays. Where true, this is the region
                             used for averaging. The whole domain will be added to the 
                             list, or will be the only domain if left to be default.
     region_names (list)   : List of strings, the names used for each region.
    '''
    
    # Open stats file and domain files
    ds_stats = xr.open_dataset(fn_stats, chunks={'profile':10000})
    dom = xr.open_dataset(fn_nemo_domain)

	# Load stats file
    ds_stats.load()
    
    print(ds_stats) 
    # Restrict time if required or define start and end dates
    if start_date is not None:
        t_ind = pd.to_datetime( ds_stats.time.values ) >= start_date
        ds_stats = ds_stats.isel(profile=t_ind)
    else:
        start_date = min(ds_stats.time)
        
    if end_date is not None:
        t_ind = pd.to_datetime( ds_stats.time.values ) <= start_date
        ds_stats = ds_stats.isel(profile=t_ind)
    else:
        start_date = min(ds_stats.time)
    
    bathymetry=False 
    # Were bottom errors calculated?
    if 'sbt_err' in ds_stats.keys():
        bathymetry=True
    
    # Some old version of files contain season dimensions, which caused problems.
    if 'season' in ds_stats.dims:
        ds_stats = ds_stats.drop_dims('season')
        
    # Make a copy of region lists to avoid infinitely growing lists..
    regional_masks = regional_masks.copy()
    region_names = region_names.copy()
    
    # Add whole domain to regions
    regional_masks.append(np.ones(dom.glamt.values.squeeze().shape))
    region_names.append('whole_domain')
    
    # Calculate the nearest indices - DONE IN ANALYSIS ROUTINE NO LONGER NEEDED
    n_regions = len(regional_masks)
    #ind2D = gu.nearest_indices_2D(dom.glamt.values.squeeze(), dom.gphit.values.squeeze(),
    #                              ds_stats.longitude.values, ds_stats.latitude.values)
    # Determine which region each point lies in
    is_in_region = [mm[ds_stats.nn_ind_y, ds_stats.nn_ind_x] for mm in regional_masks]
    is_in_region = np.array(is_in_region, dtype=bool)
    
    # 5 Seasons, define array
    reg_array = np.zeros((n_regions, 5))*np.nan
    
    # Define array to contain averaged data
    ds_mean = xr.Dataset(coords = dict(
                        longitude = ("profile", ds_stats.longitude.values),
                        latitude = ("profile", ds_stats.latitude.values),
                        time = ("profile", ds_stats.time.values),
                        region_names = ('region', region_names)),
                    data_vars = dict(
                        sst_me = (["region", "season"],  reg_array.copy()),
                        sss_me = (["region", "season"],  reg_array.copy()),
                        sst_mae = (["region", "season"], reg_array.copy()),
                        sss_mae = (["region", "season"], reg_array.copy()),
                        sst_estd = (["region", "season"],  reg_array.copy()),
                        sss_estd = (["region", "season"],  reg_array.copy()),
                        sst_crps2_mean = (["region", "season"], reg_array.copy()),
                        sss_crps2_mean = (["region", "season"], reg_array.copy()),
                        sst_crps4_mean = (["region", "season"], reg_array.copy()),
                        sss_crps4_mean = (["region", "season"], reg_array.copy()),
                        sst_crps6_mean = (["region", "season"], reg_array.copy()),
                        sss_crps6_mean = (["region", "season"], reg_array.copy())))
    
    if bathymetry:
        ds_mean['sbt_me'] = (['region','season'], reg_array.copy())
        ds_mean['sbs_me'] = (['region','season'], reg_array.copy())
        ds_mean['sbt_estd'] = (['region','season'], reg_array.copy())
        ds_mean['sbs_estd'] = (['region','season'], reg_array.copy())
        ds_mean['sbt_mae'] = (['region','season'], reg_array.copy())
        ds_mean['sbs_mae'] = (['region','season'], reg_array.copy())                                          
    
    season_indices = {'DJF':0, 'JJA':1, 'MAM':2, 'SON':3}
 
    # Loop over regions. For each, group into seasons and average.
    # Place into ds_mean dataset.
    
    bad_flag = ds_stats.bad_flag.values
    ds_stats_clean = ds_stats.isel(profile = bad_flag == False)
    is_in_region_clean = is_in_region[:, bad_flag == False]
    
    for reg in range(0,n_regions):
        reg_ind = np.where( is_in_region_clean[reg].astype(bool) )[0]
        
        if len(reg_ind)<1:
            continue

        ds_reg = ds_stats_clean.isel(profile = reg_ind)
        ds_reg_group = ds_reg.groupby('time.season')
        
        # MEANS
        ds_reg_mean = ds_reg_group.mean(dim = 'profile', skipna=True).compute()
        
        s_in_mean = ds_reg_mean.season.values
        s_ind = np.array([season_indices[ss] for ss in s_in_mean], dtype=int)

        print(s_ind)    
        ds_mean['sst_me'][reg, s_ind]  = ds_reg_mean.sst_err.values
        ds_mean['sss_me'][reg, s_ind]  = ds_reg_mean.sss_err.values
        ds_mean['sst_mae'][reg, s_ind] = ds_reg_mean.sst_abs_err.values
        ds_mean['sss_mae'][reg, s_ind] = ds_reg_mean.sss_abs_err.values
        ds_mean['sst_crps2_mean'][reg, s_ind] = ds_reg_mean.sst_crps2.values
        ds_mean['sss_crps2_mean'][reg, s_ind] = ds_reg_mean.sss_crps2.values
        ds_mean['sst_crps4_mean'][reg, s_ind] = ds_reg_mean.sst_crps4.values
        ds_mean['sss_crps4_mean'][reg, s_ind] = ds_reg_mean.sss_crps4.values
        ds_mean['sst_crps6_mean'][reg, s_ind] = ds_reg_mean.sst_crps6.values
        ds_mean['sss_crps6_mean'][reg, s_ind] = ds_reg_mean.sss_crps6.values
        
        if bathymetry:
            ds_mean['sbt_me'][reg, s_ind]  = ds_reg_mean.sbt_err.values
            ds_mean['sbs_me'][reg, s_ind]  = ds_reg_mean.sbs_err.values
            ds_mean['sbt_mae'][reg, s_ind] = ds_reg_mean.sbt_abs_err.values
            ds_mean['sbs_mae'][reg, s_ind] = ds_reg_mean.sbs_abs_err.values
        
        ds_reg_mean = ds_reg.mean(dim='profile', skipna=True).compute()
        ds_mean['sst_me'][reg, 4]  = ds_reg_mean.sst_err.values
        ds_mean['sss_me'][reg, 4]  = ds_reg_mean.sss_err.values
        ds_mean['sst_mae'][reg, 4] = ds_reg_mean.sst_abs_err.values
        ds_mean['sss_mae'][reg, 4] = ds_reg_mean.sss_abs_err.values
        ds_mean['sst_crps2_mean'][reg, 4] = ds_reg_mean.sst_crps2.values
        ds_mean['sss_crps2_mean'][reg, 4] = ds_reg_mean.sss_crps2.values
        ds_mean['sst_crps4_mean'][reg, 4] = ds_reg_mean.sst_crps4.values
        ds_mean['sss_crps4_mean'][reg, 4] = ds_reg_mean.sss_crps4.values
        ds_mean['sst_crps6_mean'][reg, 4] = ds_reg_mean.sst_crps6.values
        ds_mean['sss_crps6_mean'][reg, 4] = ds_reg_mean.sss_crps6.values
        
        if bathymetry:
            ds_mean['sbt_me'][reg, 4]  = ds_reg_mean.sbt_err.values
            ds_mean['sbs_me'][reg, 4]  = ds_reg_mean.sbs_err.values
            ds_mean['sbt_mae'][reg, 4] = ds_reg_mean.sbt_abs_err.values
            ds_mean['sbs_mae'][reg, 4] = ds_reg_mean.sbs_abs_err.values


        # STD DEVIATIONS
        ds_reg_std = ds_reg_group.std(dim = 'profile', skipna=True).compute()

        s_in_std = ds_reg_std.season.values
        s_ind = np.array([season_indices[ss] for ss in s_in_std], dtype=int)
        
        ds_mean['sst_estd'][reg, s_ind]  = ds_reg_std.sst_err.values
        ds_mean['sss_estd'][reg, s_ind]  = ds_reg_std.sss_err.values

        if bathymetry:
            ds_mean['sbt_estd'][reg, s_ind]  = ds_reg_std.sbt_err.values
            ds_mean['sbs_estd'][reg, s_ind]  = ds_reg_std.sbs_err.values

        ds_reg_std = ds_reg.std(dim='profile', skipna=True).compute()
        ds_mean['sst_estd'][reg, 4]  = ds_reg_std.sst_err.values
        ds_mean['sss_estd'][reg, 4]  = ds_reg_std.sss_err.values

        if bathymetry:
            ds_mean['sbt_estd'][reg, 4]  = ds_reg_std.sbt_err.values
            ds_mean['sbs_estd'][reg, 4]  = ds_reg_std.sbs_err.values
    
    ds_mean['start_date'] = start_date
    ds_mean['end_date'] = end_date
    ds_mean['is_in_region'] = (['region', 'profile'], is_in_region)
    ds_mean['bad_flag'] = (['profile'], ds_stats.bad_flag.values)
    
    # Write to file    
    write_ds_to_file(ds_mean, fn_out)

class analyse_and_extract():
    
    def __init__(self, fn_nemo_data, fn_nemo_domain, fn_en4, fn_out, 
                 surface_def=2, bottom_def=10,
                 nemo_chunks={'time_counter':50},
                 bathymetry = None, dist_crit=5,
                 start_date = None, end_date = None):
                 
        ''' 
        Analysis of hourly model temperature and salinity output with EN4 profiles.
        
        INPUT
         fn_nemo_data (str)   : Absolute path to output files containing hourly data
         fn_nemo_domain (str) : Absolute path to nemo domain file
         fn_en4 (str)         : Absolute path to EN4 profile files
         fn_out (str)         : Absolute path to desired output file (1 file)
         surface_def (float)  : Definition of the 'surface' in metres - for averaging
         bottom_def (float)   : Definition of the 'bottom' in metres - for averaging
         bathymetry (2Darray) : Bathymetry data to use for bottom definition.
                                If not supplied, only surface will be analysed
        '''
        
        # Open NEMO data files and define some arrays
        nemo = coast.NEMO(fn_nemo_data, fn_nemo_domain, multiple=True, chunks=nemo_chunks)
        if 'bottom_level' in nemo.dataset:
            nemo_mask = nemo.dataset.bottom_level == 0
        else:
            dom = xr.open_dataset(fn_nemo_domain)
            nemo_mask = dom.mbathy.squeeze() == 0
        nemo.dataset = nemo.dataset.rename({'t_dim':'time'})
        if bathymetry is not None:
            nemo.dataset = nemo.dataset[['votemper_top','vosaline_top',
                                         'votemper_bot','vosaline_bot']]
        else:
            nemo.dataset = nemo.dataset[['votemper_top','vosaline_top']]
        
        # Open EN4 data files
        en4 = coast.PROFILE()
        en4.read_EN4(fn_en4, multiple=True)
        
        # Cut out just the EN4 profiles inside model domain
        lonmax = np.nanmax(nemo.dataset['longitude'])
        lonmin = np.nanmin(nemo.dataset['longitude'])
        latmax = np.nanmax(nemo.dataset['latitude'])
        latmin = np.nanmin(nemo.dataset['latitude'])
        ind = coast.general_utils.subset_indices_lonlat_box(en4.dataset['longitude'], 
                                                            en4.dataset['latitude'],
                                                            lonmin, lonmax, 
                                                            latmin, latmax)[0]
        en4 = en4.isel(profile=ind)
        
        # Restrict time if required or define start and end dates
        nemo.dataset.time.load()
        if start_date is not None:
            t_ind = pd.to_datetime( nemo.dataset.time.values ) >= start_date
            nemo.dataset = nemo.dataset.isel(profile=t_ind)
        else:
            start_date = min(nemo.dataset.time)
            
        if end_date is not None:
            t_ind = pd.to_datetime( nemo.dataset.time.values ) <= start_date
            nemo.dataset = nemo.dataset.isel(profile=t_ind)
        else:
            start_date = min(nemo.dataset.time)
        
        # Cut out obs inside model time window
        n_nemo_time = nemo.dataset.dims['time']
        en4.dataset.time.load()
        nemo_time = pd.to_datetime(nemo.dataset.time.values)
        en4_time = pd.to_datetime(en4.dataset.time.values)
        time_max = max(nemo_time) + timedelta(hours=1)
        time_min = min(nemo_time) - timedelta(hours=1)
        time_ind0 = en4_time <= time_max
        time_ind1 = en4_time >= time_min
        time_ind = np.logical_and(time_ind0, time_ind1)
        en4 = en4.isel(profile=time_ind)
        
        # Get nearest model indices to observations
        en4_time = pd.to_datetime(en4.dataset.time.values)
        ind2D = gu.nearest_indices_2D(nemo.dataset.longitude, nemo.dataset.latitude, 
                                      en4.dataset.longitude, en4.dataset.latitude,
                                      mask=nemo_mask)
        
        # Vector of bools which save whether a datapoint is bad for some reason
        n_prof = len(en4_time)
        bad_flag = np.zeros(n_prof).astype(bool)
        
        # Determine if nearest points are further than dist crit away
        # If so, add bad flag
        mod_lon = nemo.dataset.longitude.isel(x_dim = ind2D[0], y_dim=ind2D[1]).values
        mod_lat = nemo.dataset.latitude.isel(x_dim=ind2D[0], y_dim=ind2D[1]).values
        obs_lon = en4.dataset.longitude.values
        obs_lat = en4.dataset.latitude.values
        interp_dist = gu.calculate_haversine_distance( mod_lon, mod_lat, 
                                                       obs_lon, obs_lat )
        bad_ind = interp_dist > dist_crit
        bad_flag[bad_ind] = True
        
        # Estimate EN4 SST as mean of top levels
        surface_ind = en4.dataset.depth <= surface_def
        
        sst_en4 = en4.dataset.potential_temperature.where(surface_ind, np.nan)
        sss_en4 = en4.dataset.practical_salinity.where(surface_ind, np.nan)
        
        sst_en4 = sst_en4.mean(dim="z_dim", skipna=True).load()
        sss_en4 = sss_en4.mean(dim="z_dim", skipna=True).load()
        
        # Bottom values
        if bathymetry is not None:
            bathy_pts = bathymetry.isel(x_dim = ind2D[0], y_dim = ind2D[1]).swap_dims({'dim_0':'profile'})
            bottom_ind = en4.dataset.depth >= (bathy_pts - bottom_def)

            sbt_en4 = en4.dataset.potential_temperature.where(bottom_ind, np.nan)
            sbs_en4 = en4.dataset.practical_salinity.where(bottom_ind, np.nan)
        
            sbt_en4 = sbt_en4.mean(dim="z_dim", skipna=True).load()
            sbs_en4 = sbs_en4.mean(dim="z_dim", skipna=True).load()
        
        # For every EN4 profile, determine the nearest model time index
        # If more than t_crit away from nearest, then discard it
        
        # Define analysis arrays
        n_prof = en4.dataset.dims['profile']
        
        sst_e = np.zeros(n_prof)*np.nan
        sss_e = np.zeros(n_prof)*np.nan
        sst_ae = np.zeros(n_prof)*np.nan
        sss_ae = np.zeros(n_prof)*np.nan
        crps_tem_2 = np.zeros(n_prof)*np.nan
        crps_sal_2 = np.zeros(n_prof)*np.nan
        crps_tem_4 = np.zeros(n_prof)*np.nan
        crps_sal_4 = np.zeros(n_prof)*np.nan
        crps_tem_6 = np.zeros(n_prof)*np.nan
        crps_sal_6 = np.zeros(n_prof)*np.nan
        
        sbt_e = np.zeros(n_prof)*np.nan
        sbs_e = np.zeros(n_prof)*np.nan
        sbt_e = np.zeros(n_prof)*np.nan
        sbs_e = np.zeros(n_prof)*np.nan
        
        x_dim_len = nemo.dataset.dims['x_dim']
        y_dim_len = nemo.dataset.dims['y_dim']
        
        n_r = nemo.dataset.dims['y_dim']
        n_c = nemo.dataset.dims['x_dim']
        n_season = 5
        
        print('Starting analysis')
        
        # LOOP over model time snapshots -- For each snapshot identify all 
        # corresponding EN4 datapoints
        for tii in range(0, n_nemo_time):
            
            print(nemo_time[tii], flush=True)
            
            
            time_diff = np.abs( nemo_time[tii] - en4_time ).astype('timedelta64[m]')
            use_ind = np.where( time_diff.astype(int) < 30 )[0]
            n_use = len(use_ind)
            
            if n_use>0:
                
                # Index the model data to time slice and load it.
                tmp = nemo.isel(time = tii).dataset
                tmp.load()
                x_tmp = ind2D[0][use_ind]
                y_tmp = ind2D[1][use_ind]
                
                x_tmp = xr.where(x_tmp<x_dim_len-7, x_tmp, np.nan)
                y_tmp = xr.where(y_tmp<y_dim_len-7, y_tmp, np.nan)

                x_tmp = xr.where(x_tmp>7, x_tmp, np.nan)
                y_tmp = xr.where(y_tmp>7, y_tmp, np.nan)
                
                shared_mask = np.logical_or(np.isnan(x_tmp), np.isnan(y_tmp))
                shared_mask = np.where(~shared_mask)
                
                x_tmp = x_tmp[shared_mask].astype(int)
                y_tmp = y_tmp[shared_mask].astype(int)
                use_ind = use_ind[shared_mask].astype(int)
                
                n_use = len(use_ind)
                if n_use<1:
                    continue
                
                # Surface errors
                tmp_pts = tmp.isel(x_dim = x_tmp, y_dim = y_tmp)
                sst_en4_tmp = sst_en4.values[use_ind]
                sss_en4_tmp = sss_en4.values[use_ind]
                sst_e[use_ind] = tmp_pts.votemper_top.values - sst_en4_tmp
                sss_e[use_ind] = tmp_pts.vosaline_top.values - sss_en4_tmp
                
                # Bottom errors
                if bathymetry is not None:
                    sbt_en4_tmp = sbt_en4.values[use_ind]
                    sbs_en4_tmp = sbs_en4.values[use_ind]
                    sbt_e[use_ind] = tmp_pts.votemper_bot.values - sbt_en4_tmp
                    sbs_e[use_ind] = tmp_pts.vosaline_bot.values - sbs_en4_tmp
                
                # CRPS 2 point radius
                nh_x = [np.arange( x_tmp[ii]-2, x_tmp[ii]+3 ) for ii in range(0,n_use)] 
                nh_y = [np.arange( y_tmp[ii]-2, y_tmp[ii]+3 ) for ii in range(0,n_use)]   
                nh = [tmp.isel(x_dim = nh_x[ii], y_dim = nh_y[ii]) for ii in range(0,n_use)] 
                crps_tem_tmp = [ cu.crps_empirical(nh[ii].votemper_top.values.flatten(), sst_en4_tmp[ii]) for ii in range(0,n_use)]
                crps_sal_tmp = [ cu.crps_empirical(nh[ii].vosaline_top.values.flatten(), sss_en4_tmp[ii]) for ii in range(0,n_use)]
                crps_tem_2[use_ind] = crps_tem_tmp
                crps_sal_2[use_ind] = crps_sal_tmp
                
                # CRPS 4 points radius
                nh_x = [np.arange( x_tmp[ii]-4, x_tmp[ii]+5 ) for ii in range(0,n_use)] 
                nh_y = [np.arange( y_tmp[ii]-4, y_tmp[ii]+5 ) for ii in range(0,n_use)]   
                nh = [tmp.isel(x_dim = nh_x[ii], y_dim = nh_y[ii]) for ii in range(0,n_use)] 
                crps_tem_tmp = [ cu.crps_empirical(nh[ii].votemper_top.values.flatten(), sst_en4_tmp[ii]) for ii in range(0,n_use)]
                crps_sal_tmp = [ cu.crps_empirical(nh[ii].vosaline_top.values.flatten(), sss_en4_tmp[ii]) for ii in range(0,n_use)]
                crps_tem_4[use_ind] = crps_tem_tmp
                crps_sal_4[use_ind] = crps_sal_tmp
                
                # CRPS 6 points radius
                nh_x = [np.arange( x_tmp[ii]-6, x_tmp[ii]+7 ) for ii in range(0,n_use)] 
                nh_y = [np.arange( y_tmp[ii]-6, y_tmp[ii]+7 ) for ii in range(0,n_use)]   
                nh = [tmp.isel(x_dim = nh_x[ii], y_dim = nh_y[ii]) for ii in range(0,n_use)] 
                crps_tem_tmp = [ cu.crps_empirical(nh[ii].votemper_top.values.flatten(), sst_en4_tmp[ii]) for ii in range(0,n_use)]
                crps_sal_tmp = [ cu.crps_empirical(nh[ii].vosaline_top.values.flatten(), sss_en4_tmp[ii]) for ii in range(0,n_use)]
                crps_tem_6[use_ind] = crps_tem_tmp
                crps_sal_6[use_ind] = crps_sal_tmp
                    
        print('Profile analysis done', flush=True)
        
        # Define absolute errors
        sst_ae = np.abs(sst_e)
        sss_ae = np.abs(sss_e)
        sbt_ae = np.abs(sbt_e)
        sbs_ae = np.abs(sbs_e)
        # Put everything into xarray dataset
        en4_season = get_season_index(sst_en4.time.values)
        
        ds = xr.Dataset(coords = dict(
                            longitude = ("profile", sst_en4.longitude.values),
                            latitude = ("profile", sst_en4.latitude.values),
                            time = ("profile", sst_en4.time.values),
                            season_ind = ("profile", en4_season)),
                        data_vars = dict(
                            obs_sst = ('profile', sst_en4.values),
                            obs_sss = ('profile', sss_en4.values),
                            sst_err = ("profile", sst_e),
                            sss_err = ("profile", sss_e),
                            sst_abs_err = ("profile", sst_ae),
                            sss_abs_err = ("profile", sss_ae),
                            sst_crps2 = ("profile", crps_tem_2),
                            sss_crps2 = ("profile", crps_sal_2),
                            sst_crps4 = ("profile", crps_tem_4),
                            sss_crps4 = ("profile", crps_sal_4),
                            sst_crps6 = ("profile", crps_tem_6),
                            sss_crps6 = ("profile", crps_sal_6),
                            nn_ind_x = ("profile", ind2D[0]),
                            nn_ind_y = ("profile", ind2D[1])))
        
        season_names = ['All','DJF','MAM','JJA','SON']
        ds = ds.chunk({'profile':10000})
        
        if bathymetry is not None:                        
            ds['obs_sbt'] = (['profile'], sbt_en4.values)
            ds['obs_sbs'] = (['profile'], sbs_en4.values)
            ds['sbt_err'] = (['profile'], sbt_e)
            ds['sbs_err'] = (['profile'], sbs_e)
            ds['sbt_abs_err'] = (['profile'], sbt_ae)
            ds['sbs_abs_err'] = (['profile'], sbs_ae)                                 
            
        ds_out = ds
        ds_out['bad_flag'] = (['profile'], bad_flag)
        
        # Write to file
        write_ds_to_file(ds_out, fn_out)
        
class radius_means():
    def __init__(self, fn_stats, fn_out, grid_lon, grid_lat, 
                 radius=10):
                 
        '''
        For a set of specified longitudes and latitudes, will calculate the mean of all
        statistics within a specified radius. This will have a smoothing effect on 
        profile data (horizontal smoothing).
        
        INPUTS
         fn_stats (str)  : Absolute path to statistics/extracted data file
         fn_out (str)    : Absolute path to desired output file
         grid_lon (array): 1D array containing longitude values of grid
         grid_lat (array): 1D array contain latitude values of grid
         radius (float)  : Radius over which to average in km
        '''
        
        # Open stats file
        stats_profile = xr.open_mfdataset(fn_stats, chunks={'profile':10000})
        
        # Define some seasonal lists
        seasons = ['Annual','DJF','MAM','JJA','SON']
        n_seasons = len(seasons)
        
        # Define grid and stack
        lon2, lat2 = np.meshgrid(grid_lon, grid_lat)
        lon = lon2.flatten()
        lat = lat2.flatten()
        
        # Define output arrays
        sst_err = np.zeros((n_seasons,len(lon)))
        sss_err = np.zeros((n_seasons,len(lon)))
        sst_N = np.zeros((n_seasons,len(lon)))
        sss_N = np.zeros((n_seasons,len(lon)))
        sbt_err = np.zeros((n_seasons,len(lon)))
        sbs_err = np.zeros((n_seasons,len(lon)))
        sbt_N = np.zeros((n_seasons,len(lon)))
        sbs_N = np.zeros((n_seasons,len(lon)))
        
        # Loop over each season and average into radii
        for season in range(1,5):
            ind_season = stats_profile.season_ind == season
            tmp = stats_profile.isel(profile=ind_season)
            tmp = tmp[['sst_err', 'sss_err','sbt_err','sbs_err']]
            tmp.load()
            
            # Remove outliers
            std = tmp.std(skipna = True)
            tmp['sst_err'] = xr.where(uf.fabs(tmp['sst_err']) > 5*std.sst_err, np.nan, tmp['sst_err'] )
            tmp['sbt_err'] = xr.where(uf.fabs(tmp['sbt_err']) > 5*std.sbt_err, np.nan, tmp['sbt_err'] )
            tmp['sss_err'] = xr.where(uf.fabs(tmp['sss_err']) > 5*std.sss_err, np.nan, tmp['sss_err'] )
            tmp['sbs_err'] = xr.where(uf.fabs(tmp['sbs_err']) > 5*std.sbs_err, np.nan, tmp['sbs_err'] )
            
            ind = gu.subset_indices_by_distance_BT(tmp.longitude, tmp.latitude, 
                                                lon, lat, radius=radius)
            
            # Surface variables
            tem = [tmp.sst_err.isel(profile=ii).values for ii in ind]
            sal = [tmp.sss_err.isel(profile=ii).values for ii in ind]
            sst_err[season] = [np.nanmean(temii) for temii in tem]
            sss_err[season] = [np.nanmean(salii) for salii in sal]
            sst_N[season] = [np.sum( ~np.isnan(temii) ) for temii in tem]
            sss_N[season] = [np.sum( ~np.isnan(salii) ) for salii in sal]
            
            # Bottom variables
            tem = [tmp.sbt_err.isel(profile=ii).values for ii in ind]
            sal = [tmp.sbs_err.isel(profile=ii).values for ii in ind]
            sbt_err[season] = [np.nanmean(temii) for temii in tem]
            sbs_err[season] = [np.nanmean(salii) for salii in sal]
            sbt_N[season] = [np.sum( ~np.isnan(temii) ) for temii in tem]
            sbs_N[season] = [np.sum( ~np.isnan(salii) ) for salii in sal]
                
        # Place analysis into dataset and write to output file
        ds = xr.Dataset(coords = dict(
                            longitude = ('location',lon),
                            latitude = ('location', lat),
                            season = ('season', seasons)),
                        data_vars = dict(
                            sst_err = (['season','location'], sst_err),
                            sss_err = (['season','location'], sss_err),
                            sst_N = (['season','location'], sst_N),
                            sss_N = (['season','location'], sss_N),
                            sbt_err = (['season','location'], sbt_err),
                            sbs_err = (['season','location'], sbs_err),
                            sbt_N = (['season','location'], sbt_N),
                            sbs_N = (['season','location'], sbs_N)))
        ds.to_netcdf(fn_out)

