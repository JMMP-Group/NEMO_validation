"""
@author: David Byrne (dbyrne@noc.ac.uk)
v1.2 (10-03-2021)

A set of routines for comparing daily mean model data with EN4 profile data. 
The input NEMO data should be daily mean data in MONTHLY files, with all
depths, and monthly EN4 data as downloaded from:

https://www.metoffice.gov.uk/hadobs/en4/

This script will read in one month at a time. Filenames for the NEMO and EN4
data are generated using the make_nemo_filename() and make_en4_filename()
routines. It is worth checking these functions to make sure that they adhere
to your filenames.

See the Github Wiki for more information:

https://github.com/JMMP-Group/NEMO_validation
        
"""
# UNCOMMENT IF USING A DEVELOPMENT VERSION OF COAST
import sys
sys.path.append('/home/users/dbyrne/code/COAsT/')

import coast
import coast.general_utils as coastgu
import coast.plot_util as pu
import numpy as np
import pandas as pd
import gsw
import xarray as xr
import sys
import os
import os.path
import glob
import scipy.stats as spst
import xarray.ufuncs as uf
import matplotlib.pyplot as plt
from dask.diagnostics import ProgressBar
from dateutil.relativedelta import relativedelta

##########################################################################################
### Analysis ROUTINES
##########################################################################################
 
def write_ds_to_file(ds, fn, **kwargs):
    ''' 
    Simple netcdf writing routine which checks if file already exists first 
    '''
    if os.path.exists(fn):
        os.remove(fn)
    ds.to_netcdf(fn, **kwargs)

def analyse_ts_regional(fn_nemo_domain, fn_extracted, fn_out, ref_depth,
                        ref_depth_method = 'bin',
                        regional_masks=[], region_names=[],
                        start_date = None, end_date = None):
    '''
    Routine for doing REGIONAL averaging of analysis files outputted using 
    analyse_ts_per_file(). INPUT is the output from the analysis and OUTPUT is a file 
    containing REGIONAL averaged statistics.

    INPUTS:
     fn_nemo_domain. (str)   : Absolute path to NEMO domain_cfg.nc
     fn_extracted    (str)   : Absolute path to single analysis file
     fn_out          (str)   : Absolute path to desired output file
     ref_depth.      (array) : 1D array describing the reference depths onto which model
                               and observations will be interpolated
     regional_masks. (list)  : List of 2D bool arrays describing averaging regions. Each 
                               array should have same shape as model domain. True 
                               indicates a model point within the defined region. Will 
                               always do a 'Whole Domain' region, even when undefined.
     region_names.   (list)  : List of strings. Names for each region in regional masks.
                            
    OUTPUTS
     Writes averages statistics to file.
     
    '''
    # Cast inputs to numpy array
    ref_depth = np.array(ref_depth)
    
    # Open dataset containing extracted profiles and statistics
    ds_ext = xr.open_mfdataset(fn_extracted, chunks={'profile':10000})
    
    # Open domain file and get bathymetry.
    dom = xr.open_dataset(fn_nemo_domain)
    bath = dom.hbatt.values.squeeze()
    dom.close()
    
    # Make a copy of regional masks to avoid any memory leaking type loops
    regional_masks = regional_masks.copy()
    region_names = region_names.copy()
    
    # Append whole domain mask
    regional_masks.append(np.ones(bath.shape))
    region_names.append('whole_domain')
    
    # Get numbers for array sizes
    n_regions = len(regional_masks)
    n_profiles = ds_ext.dims['profile']
    
    # Load only those variables that we want for interpolating to reference depths.
    ds = ds_ext[['mod_tem','obs_tem','mod_sal','obs_sal','obs_z', 'nn_ind_x', 'nn_ind_y', 'bad_flag']].astype('float32')
    ds.load()
    
    # Restrict time if required or define start and end dates
    if start_date is not None:
        t_ind = pd.to_datetime( ds.time.values ) >= start_date
        ds = ds.isel(profile=t_ind)
    else:
        start_date = min(ds.time)
        
    if end_date is not None:
        t_ind = pd.to_datetime( ds.time.values ) <= end_date
        ds = ds.isel(profile=t_ind)
    else:
        end_date = min(ds.time)
   
    # Update number of profiles
    n_profiles = ds.dims['profile']
    bathy_pts = bath[ds.nn_ind_y.values.astype(int), ds.nn_ind_x.values.astype(int)]
    is_in_region = [mm[ds.nn_ind_y.values.astype(int), ds.nn_ind_x.values.astype(int)] for mm in regional_masks]
    is_in_region = np.array(is_in_region, dtype=bool)
    
    #Figure out ref depths or bins to create interp array
    if ref_depth_method == 'interp':
        n_ref_depth = len(ref_depth)
    if ref_depth_method == 'bin':
        n_ref_depth = len(ref_depth) - 1
        bin_widths = np.array( [ref_depth[ii+1] - ref_depth[ii] for ii in np.arange(0,n_ref_depth)] )
        bin_mids = np.array( [ref_depth[ii] + .5*ref_depth[ii] for ii in np.arange(0,n_ref_depth)] )
    
    # Create dataset for interpolation
    ds_interp = xr.Dataset(coords = dict(
                               ref_depth = ('ref_depth', ref_depth),
                               longitude = ('profile', ds.longitude.values),
                               latitude = ('profile', ds.latitude.values),
                               time = ('profile', ds.time.values),
                               region = ('region', region_names)),
                           data_vars = dict(
                               bathy = ('profile', bathy_pts),
                               mod_tem = (['profile','ref_depth'], np.zeros((n_profiles, n_ref_depth), dtype='float32')*np.nan),
                               mod_sal = (['profile','ref_depth'], np.zeros((n_profiles, n_ref_depth), dtype='float32')*np.nan),
                               obs_tem = (['profile','ref_depth'], np.zeros((n_profiles, n_ref_depth), dtype='float32')*np.nan),
                               obs_sal = (['profile','ref_depth'], np.zeros((n_profiles, n_ref_depth), dtype='float32')*np.nan)))
    
    # INTERP1 = interpolate the obs and model to reference depths from the
    # OBS depths. Model already interpolated in extract routine.
    if ref_depth_method=='interp':
        for pp in range(0, n_profiles):
            prof = ds.isel(profile=pp).swap_dims({'level':'obs_z'}).dropna(dim='obs_z')
            if prof.dims['obs_z']>1:
                try:
                    print(pp)
                    prof_interp = prof.interp(obs_z = ref_depth)
                    dep_len = prof_interp.dims['obs_z']
                    ds_interp['mod_tem'][pp, :dep_len] = prof_interp.mod_tem.values
                    ds_interp['mod_sal'][pp, :dep_len] = prof_interp.mod_sal.values
                    ds_interp['obs_tem'][pp, :dep_len] = prof_interp.obs_tem.values
                    ds_interp['obs_sal'][pp, :dep_len] = prof_interp.obs_sal.values
                except:
                    print('{0}^^'.format(pp))
                    ds_interp['bathy'][pp] = np.nan
            else:
                print('{0}**'.format(pp))
                ds_interp['bathy'][pp] = np.nan
                
    # BIN = Bin into depth bins rather than interpolate.
    elif ref_depth_method=='bin':
        for pp in range(0, n_profiles):
            prof = ds.isel(profile=pp).swap_dims({'level':'obs_z'}).dropna(dim='obs_z')
            if prof.dims['obs_z']>1:
                try:
                    print(pp)
                    prof_interp = prof.interp(obs_z = ref_depth)
                    dep_len = prof_interp.dims['obs_z']
                    ds_interp['mod_tem'][pp, :dep_len] = prof_interp.mod_tem.values
                    ds_interp['mod_sal'][pp, :dep_len] = prof_interp.mod_sal.values
                    ds_interp['obs_tem'][pp, :dep_len] = prof_interp.obs_tem.values
                    ds_interp['obs_sal'][pp, :dep_len] = prof_interp.obs_sal.values
                except:
                    print('{0}^^'.format(pp))
                    ds_interp['bathy'][pp] = np.nan
            else:
                print('{0}**'.format(pp))
                ds_interp['bathy'][pp] = np.nan
        
    # Calculate errors with depth
    ds_interp['error_tem'] = (ds_interp.mod_tem - ds_interp.obs_tem).astype('float32')
    ds_interp['error_sal'] = (ds_interp.mod_sal - ds_interp.obs_sal).astype('float32')
    ds_interp['abs_error_tem'] = np.abs( (ds_interp.mod_tem - ds_interp.obs_tem).astype('float32') )
    ds_interp['abs_error_sal'] = np.abs( (ds_interp.mod_sal - ds_interp.obs_sal).astype('float32') )
    
    # Define dataset for regional averaging
    ds_reg_prof = xr.Dataset(coords = dict(
                                           region = ('region',region_names),
                                           ref_depth = ('ref_depth', ref_depth),
                                           season = ('season', ['DJF','JJA','MAM','SON','All'])))
    ds_reg_prof['prof_mod_tem'] = (['region','season','ref_depth'], np.zeros((n_regions, 5, n_ref_depth), dtype='float32')*np.nan)
    ds_reg_prof['prof_mod_sal'] = (['region','season','ref_depth'], np.zeros((n_regions, 5, n_ref_depth), dtype='float32')*np.nan)
    ds_reg_prof['prof_obs_tem'] = (['region','season','ref_depth'], np.zeros((n_regions, 5, n_ref_depth), dtype='float32')*np.nan)
    ds_reg_prof['prof_obs_sal'] = (['region','season','ref_depth'], np.zeros((n_regions, 5, n_ref_depth), dtype='float32')*np.nan)
    ds_reg_prof['prof_error_tem'] = (['region','season','ref_depth'], np.zeros((n_regions, 5, n_ref_depth), dtype='float32')*np.nan)
    ds_reg_prof['prof_error_sal'] = (['region','season','ref_depth'], np.zeros((n_regions, 5, n_ref_depth), dtype='float32')*np.nan)
    ds_reg_prof['prof_abs_error_tem'] = (['region','season','ref_depth'], np.zeros((n_regions, 5, n_ref_depth), dtype='float32')*np.nan)
    ds_reg_prof['prof_abs_error_sal'] = (['region','season','ref_depth'], np.zeros((n_regions, 5, n_ref_depth), dtype='float32')*np.nan)
    ds_reg_prof['mean_bathy'] = (['region','season'], np.zeros((n_regions, 5))*np.nan)
    
    season_str_dict = {'DJF':0,'JJA':1,'MAM':2,'SON':3}
    
    # Remove flagged points
    bad_flag = ds.bad_flag.values
    ds_interp_clean = ds_interp.isel(profile = bad_flag == False)
    is_in_region_clean = is_in_region[:, bad_flag == False]
    
    # Loop over regional arrays. Assign mean to region and seasonal means
    for reg in range(0,n_regions):
    	# Do regional average for the correct seasons
        reg_ind = np.where( is_in_region_clean[reg].astype(bool) )[0]
        if len(reg_ind)<1:
            continue
        reg_tmp = ds_interp_clean.isel(profile = reg_ind)
        reg_tmp_group = reg_tmp.groupby('time.season')
        reg_tmp_mean = reg_tmp_group.mean(dim='profile', skipna=True).compute()
        season_str = reg_tmp_mean.season.values
        season_ind = [season_str_dict[ss] for ss in season_str]
        
        ds_reg_prof['prof_mod_tem'][reg, season_ind] = reg_tmp_mean.mod_tem
        ds_reg_prof['prof_mod_sal'][reg, season_ind] = reg_tmp_mean.mod_sal
        ds_reg_prof['prof_obs_tem'][reg, season_ind] = reg_tmp_mean.obs_tem
        ds_reg_prof['prof_obs_sal'][reg, season_ind] = reg_tmp_mean.obs_sal
        ds_reg_prof['prof_error_tem'][reg, season_ind] = reg_tmp_mean.error_tem
        ds_reg_prof['prof_error_sal'][reg, season_ind] = reg_tmp_mean.error_sal
        ds_reg_prof['prof_abs_error_tem'][reg, season_ind] = reg_tmp_mean.error_tem
        ds_reg_prof['prof_abs_error_sal'][reg, season_ind] = reg_tmp_mean.error_sal
        ds_reg_prof['mean_bathy'][reg, season_ind] = reg_tmp_mean.bathy
        
        # Do regional averaging across all seasons
        reg_tmp_mean = reg_tmp.mean(dim='profile', skipna=True).compute()
        
        ds_reg_prof['prof_mod_tem'][reg, 4] = reg_tmp_mean.mod_tem
        ds_reg_prof['prof_mod_sal'][reg, 4] = reg_tmp_mean.mod_sal
        ds_reg_prof['prof_obs_tem'][reg, 4] = reg_tmp_mean.obs_tem
        ds_reg_prof['prof_obs_sal'][reg, 4] = reg_tmp_mean.obs_sal
        ds_reg_prof['prof_error_tem'][reg, 4] = reg_tmp_mean.error_tem
        ds_reg_prof['prof_error_sal'][reg, 4] = reg_tmp_mean.error_sal
        ds_reg_prof['prof_abs_error_tem'][reg, 4] = reg_tmp_mean.error_tem
        ds_reg_prof['prof_abs_error_sal'][reg, 4] = reg_tmp_mean.error_sal
        ds_reg_prof['mean_bathy'][reg, 4] = reg_tmp_mean.bathy
    
    # Drop bathy for some reason
    ds_interp = ds_interp.drop('bathy')
    
    # Merge output dataset
    ds_interp = xr.merge((ds_interp, ds_reg_prof))
    ds_interp['is_in_region'] = (['region','profile'], is_in_region)
    ds_interp['bad_flag'] = (['profile'], ds.bad_flag.values)
    
    ds_interp['start_date'] = start_date
    ds_interp['end_date'] = end_date
    
    # Write output to file
    write_ds_to_file(ds_interp, fn_out)

def extract_ts_per_file(fn_nemo_data, fn_nemo_domain, fn_en4, fn_out,  
                        run_name = 'Undefined', surface_def=5, bottom_def=10, 
                        dist_crit=5, n_obs_levels=400, 
                        model_frequency='daily', instant_data=False):
    '''
    Extracts and does some basic analysis and identification of model data at obs
    locations, times and depths. Writes extracted data to file.
    
    INPUTS:
     fn_nemo_data (str)   : Absolute filename to a monthly nemo file containing daily mean data.
     fn_nemo_domain (str) : Absolute filepath to corresponding NEMO domain_cfg.nc
     fn_en4 (str)         : Absolute filepath to monthly EN4 profile data file.
     fn_out (str)         : Absolute filepath for desired output file.
     run_name (str)       : Name of run. [default='Undefined']
     surface_def (float)  : Depth of surface for averaging surface variables [default=5m]
     bottom_def (float)   : Distance from bottom for defining bottom variable averaging [default=10m]
     dist_crit (float)    : Distance at which to omit datapoints if the resulting interpolated
                            point is too far from obs point. [default=5km]
                            
    OUTPUTS:
     Writes extracted data to file.
    '''
    
    # 1) Read NEMO, then extract desired variables
    try:   
        nemo = coast.NEMO(fn_nemo_data, fn_nemo_domain, chunks={'time_counter':1})
        dom = xr.open_dataset(fn_nemo_domain) 
        if 'bottom_level' in nemo.dataset:
            nemo.dataset['landmask'] = (('y_dim','x_dim'), nemo.dataset.bottom_level==0)
        else:
            nemo.dataset['landmask'] = (('y_dim','x_dim'), dom.mbathy.squeeze()==0)
            nemo.dataset['bottom_level'] = (('y_dim', 'x_dim'), dom.mbathy.squeeze())
        nemo.dataset['bathymetry'] = (('y_dim', 'x_dim'), dom.hbatt[0] )
        nemo = nemo.dataset[['temperature','salinity','depth_0', 'landmask','bathymetry', 'bottom_level']]
        nemo = nemo.rename({'temperature':'tem','salinity':'sal'})
        mod_time = nemo.time.values
        
    except:
        print('       !!!Problem with NEMO Read: {0}'.format(fn_nemo_data))
        return
        
    # 2) Read EN4, then extract desired variables
    try:
        # Read relevant EN4 files
        en4 = coast.PROFILE()
        en4.read_EN4(fn_en4, chunks={})
        en4 = en4.dataset[['potential_temperature','practical_salinity','depth']]
        en4 = en4.rename({'practical_salinity':'sal', 'potential_temperature':'tem'})
    except:        
        print('       !!!Problem with EN4 Read: {0}'.format(fn_en4))
        return
    
    print('1-2) Files read: \n >>> {0} \n >>> {1}'.format(fn_nemo_data, fn_en4), flush=True)
    
    # 3) Use only observations that are within model domain
    lonmax = np.nanmax(nemo['longitude'])
    lonmin = np.nanmin(nemo['longitude'])
    latmax = np.nanmax(nemo['latitude'])
    latmin = np.nanmin(nemo['latitude'])
    ind = coast.general_utils.subset_indices_lonlat_box(en4['longitude'], 
                                                        en4['latitude'],
                                                        lonmin, lonmax, 
                                                        latmin, latmax)[0]
    en4 = en4.isel(profile=ind)
    print('3) EN4 subsetted to model domain.', flush=True)
    
    # 4) Use only observations that are within model time window.
    en4_time = en4.time.values
    time_max = pd.to_datetime( max(mod_time) ) + relativedelta(hours=12)
    time_min = pd.to_datetime( min(mod_time) ) - relativedelta(hours=12)
    ind = np.logical_and( en4_time >= time_min, en4_time <= time_max )
    en4 = en4.isel(profile=ind)
    en4.load()
    print('4) EN4 subsetted to model time period.', flush=True)
    
    # ----------------------------------------------------
    # 5) Get model indices (space and time) corresponding to observations
    # Does a basic nearest neighbour analysis in time and space.
    
    # SPATIAL indices
    ind2D = coastgu.nearest_indices_2D(nemo['longitude'], nemo['latitude'],
                                       en4['longitude'], en4['latitude'], 
                                       mask=nemo.landmask)
    
    print('Spatial Indices Calculated', flush=True)
    
    # TIME indices
    en4_time = en4.time.values
    ind_time = [ np.argmin( np.abs( mod_time - en4_time[tt] ) ) for tt in range(en4.dims['profile'])]
    min_time = [ np.min( np.abs( mod_time - en4_time[tt] ) ).astype('timedelta64[h]') for tt in range(en4.dims['profile'])]
    ind_time = xr.DataArray(ind_time)
    print('Time Indices Calculated', flush=True)
    
    # INDEX the data and load
    mod_profiles = nemo.isel(x_dim=ind2D[0], y_dim=ind2D[1], t_dim=ind_time)
    mod_profiles = mod_profiles.rename({'dim_0':'profile'})
    with ProgressBar():
        mod_profiles.load()
    print('Model indexed and loaded', flush=True)
    
    # Define variable arrays for interpolated data for monthly EN4 data
    n_mod_levels = mod_profiles.dims['z_dim']
    n_prof = en4.dims['profile']
    data = xr.Dataset(coords = dict(
                          longitude=     (["profile"], en4.longitude.values),
                          latitude=      (["profile"], en4.latitude.values),
                          time=          (["profile"], en4.time.values),
                          level=         (['level'], np.arange(0,n_obs_levels)),
                          ex_longitude = (["profile"], mod_profiles.longitude.values),
                          ex_latitude =  (["profile"], mod_profiles.longitude.values),
                          ex_time =      (["profile"], mod_profiles.time.values),
                          ex_level =     (["ex_level"], np.arange(0, n_mod_levels))),
                      data_vars = dict(
                          mod_tem = (['profile','level'],  np.zeros((n_prof , n_obs_levels))*np.nan),
                          obs_tem = (['profile','level'],  np.zeros((n_prof , n_obs_levels))*np.nan),
                          mod_sal = (['profile','level'],  np.zeros((n_prof , n_obs_levels))*np.nan),
                          obs_sal = (['profile','level'],  np.zeros((n_prof , n_obs_levels))*np.nan),
                          mod_rho = (['profile','level'],  np.zeros((n_prof , n_obs_levels))*np.nan),
                          obs_rho = (['profile','level'],  np.zeros((n_prof , n_obs_levels))*np.nan),
                          mod_s0 =  (['profile', 'level'], np.zeros((n_prof , n_obs_levels))*np.nan),
                          obs_s0 =  (['profile', 'level'], np.zeros((n_prof , n_obs_levels))*np.nan),
                          mask_tem = (['profile','level'], np.zeros((n_prof , n_obs_levels))*np.nan),
                          mask_sal = (['profile','level'], np.zeros((n_prof , n_obs_levels))*np.nan),
                          mask_rho = (['profile','level'], np.zeros((n_prof , n_obs_levels))*np.nan),
                          obs_z =    (['profile','level'],    np.zeros((n_prof , n_obs_levels))*np.nan),
                          ex_mod_tem = (["profile", "ex_level"], mod_profiles.tem.values),
                          ex_mod_sal = (["profile", "ex_level"], mod_profiles.sal.values),
                          ex_depth =   (["profile", "ex_level"], mod_profiles.depth_0.values.T),
                          nn_ind_x = (["profile"], ind2D[0]),
                          nn_ind_y = (["profile"], ind2D[1])))

	# Vector of bools which save whether a datapoint is bad for some reason
    bad_flag = np.zeros(n_prof).astype(bool)
    
    # Now loop over profiles and interpolate model onto obs.
    for prof in range(0,n_prof):
        
        
        # Select the current profile
        mod_profile = mod_profiles.isel(profile = prof)
        obs_profile = en4.isel(profile = prof)
        
        # If the nearest neighbour interpolation is bad, then skip the 
        # vertical interpolation -> keep profile as nans in monthly array
        if all(np.isnan(mod_profile.tem)):
            bad_flag[prof] = True
            continue
        
        # Check that model point is within threshold distance of obs
        # If not, skip vertical interpolation -> keep profile as nans
        interp_dist = coastgu.calculate_haversine_distance(
                                             obs_profile.longitude, 
                                             obs_profile.latitude, 
                                             mod_profile.longitude, 
                                             mod_profile.latitude)
        if interp_dist > dist_crit:
            bad_flag[prof] = True
            continue
        
        # Use bottom_level to mask dry depths
        if 'bottom_level' in mod_profile:
            bl = mod_profile.bottom_level.squeeze().values
            mod_profile = mod_profile.isel(z_dim=range(0,bl))
        
        # Interpolate model to obs depths using a linear interp
        # If interpolation fails for some or any reason, skip to next iteration
        obs_profile = obs_profile.rename({'z_dim':'depth'})
        obs_profile = obs_profile.set_coords('depth')
        mod_profile = mod_profile.rename({'z_dim':'depth_0'})
        try:
            mod_profile_int = mod_profile.interp(depth_0 = obs_profile.depth.values)
        except:
            bad_flag[prof] = True
            continue
        
        # Calculate Density
        ap_obs = gsw.p_from_z( -obs_profile.depth, obs_profile.latitude )
        ap_mod = gsw.p_from_z( -obs_profile.depth, mod_profile_int.latitude )
        # Absolute Salinity            
        sa_obs = gsw.SA_from_SP( obs_profile.sal, ap_obs, 
                                obs_profile.longitude, 
                                obs_profile.latitude )
        sa_mod = gsw.SA_from_SP( mod_profile_int.sal, ap_mod, 
                                mod_profile_int.longitude, 
                                mod_profile_int.latitude )
        # Conservative Temperature
        ct_obs = gsw.CT_from_pt( sa_obs, obs_profile.tem ) 
        ct_mod = gsw.CT_from_pt( sa_mod, mod_profile_int.tem ) 
        
        # In-situ density
        obs_rho = gsw.rho( sa_obs, ct_obs, ap_obs )
        mod_rho = gsw.rho( sa_mod, ct_mod, ap_mod ) 
        
        # Potential Density
        obs_s0 = gsw.sigma0(sa_obs, ct_obs)
        mod_s0 = gsw.sigma0(sa_mod, ct_mod)
        
        # Assign to main array
        data['mod_tem'][prof] = mod_profile_int.tem.values
        data['obs_tem'][prof] = obs_profile.tem.values
        data['mod_sal'][prof] = mod_profile_int.sal.values
        data['obs_sal'][prof] = obs_profile.sal.values
        data['mod_rho'][prof] = mod_rho
        data['obs_rho'][prof] = obs_rho
        data['mod_s0'][prof] = mod_s0
        data['obs_s0'][prof] = obs_s0
        data['obs_z'][prof] = obs_profile.depth
        
    print(' Interpolated Profiles.', flush=True)
    
    # Define seasons as month numbers and identify datapoint seasons
    month_season_dict = {1:1, 2:1, 3:2, 4:2, 5:2, 6:3,
                         7:3, 8:3, 9:4, 10:4, 11:4, 12:1}
    pd_time = pd.to_datetime(data.time.values)
    pd_month = pd_time.month
    season_save = [month_season_dict[ii] for ii in pd_month]
    
    data['season'] = ('profile', season_save)
    data.attrs['run_name'] = run_name
    data['bad_flag'] = ('profile', bad_flag)
    
    # Errors at all depths
    data["error_tem"] = (['profile','level'], data.mod_tem - data.obs_tem)
    data["error_sal"] = (['profile','level'], data.mod_sal - data.obs_sal)
    
    # Absolute errors at all depths
    data["abs_error_tem"] = (['profile','level'], np.abs(data.error_tem))
    data["abs_error_sal"] = (['profile','level'], np.abs(data.error_sal))
    
    # Mean errors across depths
    data['me_tem'] = ('profile', np.nanmean(data.error_tem, axis=1))
    data['me_sal'] = ('profile', np.nanmean(data.error_sal, axis=1))
    
    # Mean absolute errors across depths
    data['mae_tem'] = (['profile'], np.nanmean(data.abs_error_tem, axis=1))
    data['mae_sal'] = (['profile'], np.nanmean(data.abs_error_tem, axis=1))
              
    # Write monthly stats to file
    write_ds_to_file(data, fn_out, mode='w', unlimited_dims='profile')
    
    print(' >>>>>>>  File Written: ' + fn_out, flush=True)
    
    return 

def concatenate_output_files(files, fn_out):
    
    fn_list = glob.glob(files)
    ds_list = [xr.open_dataset(ff, chunks={'profile':10000}) for ff in fn_list]
    ds_concat = xr.concat(ds_list, dim='profile')
    write_ds_to_file(ds_concat, fn_out)
    return
    
##########################################################################################
### PLOTTING ROUTINES
##########################################################################################

class plot_ts_monthly_single_cfg():
    
    def __init__(self, fn_profile_stats, dn_out, run_name, file_type='.png'):

        stats = xr.open_mfdataset(fn_profile_stats, chunks={})
        
        #Loop over seasons
        seasons = ['All','DJF','MAM','JJA','SON']
        lonmax = np.nanmax(stats.longitude)
        lonmin = np.nanmin(stats.longitude)
        latmax = np.nanmax(stats.latitude)
        latmin = np.nanmin(stats.latitude)
        lonbounds = [lonmin-1, lonmax+1]
        latbounds = [latmin-1, latmax+1]
        
        for ss in range(0, 5):
        
            if ss>0:
                ind_season = stats.season==ss
                stats_tmp = stats.isel(profile=ind_season)    
            else:
                stats_tmp = stats
        
            # Surface TEMPERATURE
            f,a = pu.create_geo_axes(lonbounds, latbounds)
            sca = a.scatter(stats_tmp.longitude, stats_tmp.latitude, c=stats_tmp.surf_error_tem, 
                vmin=-1.5, vmax=1.5, linewidths=0, zorder=100, cmap='seismic',s=1)
            f.colorbar(sca)
            a.set_title('Monthly EN4 SST Anom. (degC) | {0} | {1} - EN4'.format(seasons[ss], run_name), fontsize=9)
            fn_out = 'en4_surf_error_tem_{0}_{1}{2}'.format(seasons[ss], run_name, file_type)
            fn_out = os.path.join(dn_out, fn_out)
            f.savefig(fn_out)
            
            #Surface SALINITY
            f,a = pu.create_geo_axes(lonbounds, latbounds)
            sca = a.scatter(stats_tmp.longitude, stats_tmp.latitude, c=stats_tmp.surf_error_sal, 
                vmin=-1.5, vmax=1.5, linewidths=0, zorder=100, cmap='seismic', s=1)
            f.colorbar(sca)
            a.set_title('Monthly EN4 SSS Anom. (PSU) | {0} | {1} - EN4'.format(seasons[ss], run_name), fontsize=9)
            fn_out = 'en4_surf_error_sal_{0}_{1}{2}'.format(seasons[ss], run_name, file_type)
            fn_out = os.path.join(dn_out, fn_out)
            f.savefig(fn_out)
            
            # Bottom TEMPERATURE
            f,a = pu.create_geo_axes(lonbounds, latbounds)
            sca = a.scatter(stats_tmp.longitude, stats_tmp.latitude, c=stats_tmp.bott_error_tem, 
                vmin=-1.5, vmax=1.5, linewidths=0, zorder=100, cmap='seismic', s=1)
            f.colorbar(sca)
            a.set_title('Monthly EN4 SBT Anom. (degC) | {0} | {1} - EN4'.format(seasons[ss], run_name), fontsize=9)
            fn_out = 'en4_bott_error_tem_{0}_{1}{2}'.format(seasons[ss], run_name, file_type)
            fn_out = os.path.join(dn_out, fn_out)
            f.savefig(fn_out)
            
            # Bottom SALINITY
            f,a = pu.create_geo_axes(lonbounds, latbounds)
            sca = a.scatter(stats_tmp.longitude, stats_tmp.latitude, c=stats_tmp.bott_error_sal, 
                vmin=-1.5, vmax=1.5, linewidths=0, zorder=100, cmap='seismic', s=1)
            f.colorbar(sca)
            a.set_title('Monthly EN4 SBS Anom. (PSU) | {0} | {1} - EN4'.format(seasons[ss], run_name), fontsize=9)
            fn_out = 'en4_bott_error_sal_{0}_{1}{2}'.format(seasons[ss], run_name, file_type)
            fn_out = os.path.join(dn_out, fn_out)
            f.savefig(fn_out)
    
    
class plot_ts_monthly_multi_cfg():
    def __init__(self, fn_regional_stats, dn_out, run_name, file_type='.png'):
        
        stats_list = [xr.open_dataset(ff, chunks={}) for ff in fn_regional_stats]
    
        # For titles
        region_names = ['North Sea','Outer Shelf','Norwegian Trench','English Channel','Whole Domain']
        # For file names
        region_abbrev = ['northsea','outershelf','nortrench','engchannel','wholedomain']
        
        season_names = ['Annual','DJF','MAM','JJA','SON']
        legend = run_name
        n_regions = len(region_names)
        n_seasons = len(season_names)
        for rr in range(0,n_regions):
            for ss in range(1,n_seasons):
                tem_list = [tmp.prof_error_tem.isel( region=rr, season=ss, depth=np.arange(0,30) ) for tmp in stats_list]
                sal_list = [tmp.prof_error_sal.isel( region=rr, season=ss, depth=np.arange(0,30) ) for tmp in stats_list]
                
                title_tmp = '$\Delta T$ (degC) | {0} | {1}'.format(region_names[rr], season_names[ss])
                fn_out = 'prof_error_tem_{0}_{1}{2}'.format(season_names[ss], region_abbrev[rr], file_type)
                fn_out = os.path.join(dn_out, fn_out)
                f,a = self.plot_profile_centred(tem_list[0].depth, tem_list,
                              title = title_tmp, legend_names = legend)
                print("  >>>>>  Saving: " + fn_out)
                f.savefig(fn_out)
                plt.close()
                
                title_tmp = '$\Delta S$ (PSU) |' + region_names[rr] +' | '+season_names[ss]
                fn_out = 'prof_error_sal_{0}_{1}{2}'.format(season_names[ss], region_abbrev[rr], file_type)
                fn_out = os.path.join(dn_out, fn_out)
                f,a = self.plot_profile_centred(sal_list[0].depth, sal_list,
                             title = title_tmp,legend_names = legend)
                print("  >>>>>  Saving: " + fn_out)
                f.savefig(fn_out)
                plt.close()
    
                
    def plot_profile_centred(self, depth, variables, title="", legend_names= {} ):
    
        fig = plt.figure(figsize=(3.5,7))
        ax = plt.subplot(111)
        
        if type(variables) is not list:
            variables = [variables]
    
        xmax = 0
        for vv in variables:
            xmax = np.max([xmax, np.nanmax(np.abs(vv))])
            ax.plot(savgol_filter(vv.squeeze(),5,2), depth.squeeze())
            
        plt.xlim(-xmax-0.05*xmax, xmax+0.05*xmax)
        ymax = np.nanmax(np.abs(depth))
        plt.plot([0,0],[-1e7,1e7], linestyle='--',linewidth=1,color='k')
        plt.ylim(0,ymax)
        plt.gca().invert_yaxis()
        plt.ylabel('Depth (m)')
        plt.grid()
        plt.legend(legend_names, fontsize=10)
        plt.title(title, fontsize=8)
        return fig, ax
                
    def plot_profile(self, depth, variables, title, fn_out, legend_names= {} ):
    
        fig = plt.figure(figsize=(3.5,7))
        ax = plt.subplot(111)
        
        if type(variables) is not list:
            variables = [variables]
    
        for vv in variables:
            ax.plot(vv.squeeze(), depth.squeeze())
        plt.gca().invert_yaxis()
        plt.ylabel('Depth (m)')
        plt.grid()
        plt.legend(legend_names)
        plt.title(title, fontsize=12)
        print("  >>>>>  Saving: " + fn_out)
        plt.savefig(fn_out)
        plt.close()
        return fig, ax
