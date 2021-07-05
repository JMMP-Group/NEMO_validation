"""
A set of routines for comparing daily mean model data with EN4 profile data. 

The two main methods in this module are:
    
    1. extract_ts_per_file()
    2. analyse_ts_regional()

The first extracts model temperature and salinity profiles at the nearest 
corresponding EN4 location. The second Averages this data, along with errors
into regional and seasonal boxes. The second routine takes the output from
the first.

See docstrings of the routines for more information or the Github Wiki:
    
https://github.com/JMMP-Group/NEMO_validation

EN4 Data:

https://www.metoffice.gov.uk/hadobs/en4/
        
"""
# UNCOMMENT IF USING A DEVELOPMENT VERSION OF COAST
import sys
sys.path.append('/home/users/dbyrne/code/COAsT/')

import coast
import coast.general_utils as coastgu
import numpy as np
import pandas as pd
import xarray as xr
import sys
import os
import os.path
import glob
import scipy.stats as spst
import scipy.interpolate as interp
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
                        ref_depth_method = 'interp',
                        regional_masks=[], region_names=[],
                        start_date = None, end_date = None, dist_omit=100):
    '''
    VERSION 1.3 (30/06/2021)
    
    Routine for doing REGIONAL and SEASONAL averaging of analysis files outputted using 
    extract_ts_per_file(). INPUT is the output from the analysis and OUTPUT is a file 
    containing REGIONAL averaged statistics. Multiple files can be provided
    to this analysis but it will be quicker if concatenated beforehand.

    INPUTS:
     fn_nemo_domain. (str)   : Absolute path to NEMO domain_cfg.nc or mesh_mask.nc
                               Used only for bathymetry
     fn_extracted    (str)   : Absolute path to single analysis file
     fn_out          (str)   : Absolute path to desired output file
     ref_depth.      (array) : 1D array describing the reference depths onto which model
                               and observations will be interpolated
     ref_depth_method (str)  : 'interp' or 'bin'. If interp, then routine will
                               interpolate both model and observed values from
                               observation depths onto the common ref_depth. If
                               bin, then ref_depth will be treated as the 
                               boundaries of averaging bins. BIN CURRENTLY
                               UNIMPLEMENTED. [default = 'interp']
     regional_masks. (list)  : List of 2D bool arrays describing averaging regions. Each 
                               array should have same shape as model domain. True 
                               indicates a model point within the defined region. Will 
                               always do a 'Whole Domain' region, even when undefined.
     region_names.   (list)  : List of strings. Names for each region in regional masks.
                               To be saved in output file.
     start_date (datetime)   : Start date for analysis
     end_date (datetime)     : End date for analysis
     dist_omit (float)       : Distance of nearest grid cell from observation
                               at which to omit the datapoint from averaging (km)
                            
    OUTPUTS
     Writes averages statistics to file. netCDF file has dimensions:
         profile   : Profile location dimension
         ref_depth : Reference depth dimension
                     If ref_depth_method == 'bin', then this will be bin
                     midpoints.
         region    : Regional dimension
         season    : Season dimension
     Data Variables:
         mod_tem   : Model temperature on reference depths/bins
         mod_sal   : Model salinity on reference depths/bins
         obs_tem   : Observed temperature on reference depths/bins
         obs_sal   : Observed salinity on reference depths/bins
         error_tem : Temperature errors on reference depths
         error_sal : Salinity errors on reference depths
         abs_error_tem : Absolute temperature err on reference depths
         abs_error_sal : Absolute salinity err on reference depths
         prof_mod_tem  : Regional/seasonal averaged model temperatures
         prof_mod_sal  : Regional/seasonal averaged model salinity
         prof_obs_tem  : Regional/seasonal averaged observed temperature
         prof_obs_sal  : Regional/seasonal averaged obs salinity
         prof_error_tem     : Regional/seasonal averaged temperature error
         prof_error_sal     : Regional/seasonal averaged salinity error
         prof_abs_error_tem : Regional/seasonal averaged abs. temp. error
         prof_abs_error_sal : Regional/seasonal averaged abs. sal. error
         mean_bathy   : Mean model bathymetric depths for profiles used in each
                        region/season
         is_in_region : Boolean array. Described whether each profile is within
                        each region
         start_date   : Start date for analysis
         end_date     : End date for analysis
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
    ds = ds_ext[['mod_tem','obs_tem','mod_sal','obs_sal','obs_z', 'nn_ind_x', 'nn_ind_y', 'interp_dist']].astype('float32')
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
    
    # Figure out which points lie in which region
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
            print(pp)
            prof = ds.isel(profile=pp)
            
            mask = np.isnan( prof.obs_sal )
            prof['mod_sal'][mask] = np.nan
            mask = np.isnan( prof.obs_tem )
            prof['mod_tem'][mask] = np.nan
            
            depth0 = prof.obs_z.values
                
            f = interp.interp1d(depth0, prof.mod_tem.values, fill_value=np.nan, bounds_error=False)
            ds_interp['mod_tem'][pp] = f( ref_depth )
            f = interp.interp1d(depth0, prof.mod_sal.values, fill_value=np.nan, bounds_error=False)
            ds_interp['mod_sal'][pp] = f( ref_depth )
            f = interp.interp1d(depth0, prof.obs_tem.values, fill_value=np.nan, bounds_error=False)
            ds_interp['obs_tem'][pp] = f( ref_depth )
            f = interp.interp1d(depth0, prof.obs_sal.values, fill_value=np.nan, bounds_error=False)
            ds_interp['obs_sal'][pp] = f( ref_depth )
                
    # BIN = Bin into depth bins rather than interpolate - NOT USED CURRENTLY
    elif ref_depth_method=='bin':
        raise NotImplementedError()
        
    # Calculate errors with depth
    ds_interp['error_tem'] = (ds_interp.mod_tem - ds_interp.obs_tem).astype('float32')
    ds_interp['error_sal'] = (ds_interp.mod_sal - ds_interp.obs_sal).astype('float32')
    ds_interp['abs_error_tem'] = np.abs( (ds_interp.mod_tem - ds_interp.obs_tem).astype('float32') )
    ds_interp['abs_error_sal'] = np.abs( (ds_interp.mod_sal - ds_interp.obs_sal).astype('float32') )
    
    # Define dataset for regional averaging
    empty_array = np.zeros((n_regions, 5, n_ref_depth), dtype='float32')*np.nan
    ds_reg_prof = xr.Dataset(coords = dict(
                                region = ('region',region_names),
                                ref_depth = ('ref_depth', ref_depth),
                                season = ('season', ['DJF','JJA','MAM','SON','All'])),
                             data_vars = dict(
                                prof_mod_tem = (['region','season','ref_depth'], empty_array.copy()),
                                prof_mod_sal = (['region','season','ref_depth'], empty_array.copy()),
                                prof_obs_tem = (['region','season','ref_depth'], empty_array.copy()),
                                prof_obs_sal = (['region','season','ref_depth'], empty_array.copy()),
                                prof_error_tem = (['region','season','ref_depth'], empty_array.copy()),
                                prof_error_sal = (['region','season','ref_depth'], empty_array.copy()),
                                prof_abs_error_tem = (['region','season','ref_depth'], empty_array.copy()),
                                prof_abs_error_sal = (['region','season','ref_depth'], empty_array.copy()),
                                mean_bathy = (['region','season'], np.zeros((n_regions, 5))*np.nan)))
    
    season_str_dict = {'DJF':0,'JJA':1,'MAM':2,'SON':3, 'All':4}
    
    # Remove flagged points
    omit_flag = ds.interp_dist.values <= dist_omit
    ds_interp_clean = ds_interp.isel(profile = omit_flag)
    is_in_region_clean = is_in_region[:, omit_flag]
    
    # Loop over regional arrays. Assign mean to region and seasonal means
    for reg in range(0,n_regions):
        print(reg)
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
        ds_reg_prof['prof_abs_error_tem'][reg, season_ind] = reg_tmp_mean.abs_error_tem
        ds_reg_prof['prof_abs_error_sal'][reg, season_ind] = reg_tmp_mean.abs_error_sal
        ds_reg_prof['mean_bathy'][reg, season_ind] = reg_tmp_mean.bathy
        
        # Do regional averaging across all seasons
        reg_tmp_mean = reg_tmp.mean(dim='profile', skipna=True).compute()
        
        ds_reg_prof['prof_mod_tem'][reg, 4] = reg_tmp_mean.mod_tem
        ds_reg_prof['prof_mod_sal'][reg, 4] = reg_tmp_mean.mod_sal
        ds_reg_prof['prof_obs_tem'][reg, 4] = reg_tmp_mean.obs_tem
        ds_reg_prof['prof_obs_sal'][reg, 4] = reg_tmp_mean.obs_sal
        ds_reg_prof['prof_error_tem'][reg, 4] = reg_tmp_mean.error_tem
        ds_reg_prof['prof_error_sal'][reg, 4] = reg_tmp_mean.error_sal
        ds_reg_prof['prof_abs_error_tem'][reg, 4] = reg_tmp_mean.abs_error_tem
        ds_reg_prof['prof_abs_error_sal'][reg, 4] = reg_tmp_mean.abs_error_sal
        ds_reg_prof['mean_bathy'][reg, 4] = reg_tmp_mean.bathy
    
    # Drop bathy for some reason
    ds_interp = ds_interp.drop('bathy')
    
    # Merge output dataset
    ds_interp = xr.merge((ds_interp, ds_reg_prof))
    ds_interp['is_in_region'] = (['region','profile'], is_in_region)
    
    ds_interp['start_date'] = start_date
    ds_interp['end_date'] = end_date
    
    # Write output to file
    write_ds_to_file(ds_interp, fn_out)

def extract_ts_per_file(fn_nemo_data, fn_nemo_domain, fn_en4, fn_out,  
                        run_name = 'Undefined', z_interp = 'linear'):
    '''
    VERSION 1.3 (30/06/2021)
    
    Extracts and does some basic analysis and identification of model data at obs
    locations, times and depths. Writes extracted data to file. This routine
    does the analysis per file, for example monthly files. Output can be
    subsequently input to analyse_ts_regional. Data is extracted saved in two
    forms: on the original model depth levels and interpolated onto the
    observed depth levels.
    
    This routine can be used in a loop to loop over all files containing
    daily or 25hourm data. Output files can be concatenated using a tools such
    as ncks or the concatenate_output_files() routine in this module.
    
    INPUTS:
     fn_nemo_data (str)   : Absolute filename to a monthly nemo file containing 
                            daily or 25hour mean data.
     fn_nemo_domain (str) : Absolute filepath to corresponding NEMO domain_cfg.nc
     fn_en4 (str)         : Absolute filepath to monthly EN4 profile data file.
     fn_out (str)         : Absolute filepath for desired output file.
     run_name (str)       : Name of run. [default='Undefined']
                            Will be saved to output
     z_interp (str)       : Type of scipy interpolation to use for depth interpolation
                            [default = 'linear']
                            
    OUTPUTS:
     Writes extracted data to file. Extracted dataset has the dimensions:
         ex_level : Model level for directly extracted data
         level    : Observation level for interpolated model and EN4 data
         profile  : Profile location
     Output variables:
         mod_tem       : Model temperature interpolated onto obs depths
         obs_tem       : EN4 tempersture profiles
         mod_sal       : Model salinity interpolated onto obs depths
         obs_sal       : EN4 salinity profiles
         obs_z         : EN4 depths
         ex_mod_tem    : Model temperature profiles 
         ex_mod_sal    : Model salinity profiles 
         ex_depth      : Model depth at profiles
         nn_ind_x      : Nearest neighbour x (column) indices
         nn_ind_y      : Nearest neighbour y (row) indices
         season        : Season indices -> 0 = DJF, 1=JJA, 2=MAM, 3=SON
         interp_dist   : Distances from nearest model cell to EN4 locations
         error_tem     : Model-EN4 temperature differences at observation depths
         error_sal     : Model-EN4 salinity differences at observation depths
         abs_error_tem : Model-EN4 abs. temperature differences at obs depths
         abs_error_sal : Model-EN4 abs. salinity differences at obs depths
         me_tem        : Mean temperature difference over each profile
         me_sal        : Mean salinity difference over each profile
         mae_tem       : Mean absolute temperature difference over each profile
         mae_sal       : Mean absolute salinity difference over each profile
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
    n_obs_levels = en4.dims['level']
    n_prof = en4.dims['profile']
    data = xr.Dataset(coords = dict(
                          longitude=     (["profile"], en4.longitude.values),
                          latitude=      (["profile"], en4.latitude.values),
                          time=          (["profile"], en4.time.values),
                          level=         (['level'], np.arange(0,n_obs_levels)),
                          ex_longitude = (["profile"], mod_profiles.longitude.values),
                          ex_latitude =  (["profile"], mod_profiles.latitude.values),
                          ex_time =      (["profile"], mod_profiles.time.values),
                          ex_level =     (["ex_level"], np.arange(0, n_mod_levels))),
                      data_vars = dict(
                          mod_tem = (['profile','level'],  np.zeros((n_prof , n_obs_levels))*np.nan),
                          obs_tem = (['profile','level'],  np.zeros((n_prof , n_obs_levels))*np.nan),
                          mod_sal = (['profile','level'],  np.zeros((n_prof , n_obs_levels))*np.nan),
                          obs_sal = (['profile','level'],  np.zeros((n_prof , n_obs_levels))*np.nan),
                          #mod_rho = (['profile','level'],  np.zeros((n_prof , n_obs_levels))*np.nan),
                          #obs_rho = (['profile','level'],  np.zeros((n_prof , n_obs_levels))*np.nan),
                          #mod_s0 =  (['profile', 'level'], np.zeros((n_prof , n_obs_levels))*np.nan),
                          #obs_s0 =  (['profile', 'level'], np.zeros((n_prof , n_obs_levels))*np.nan),
                          #mask_tem = (['profile','level'], np.zeros((n_prof , n_obs_levels))*np.nan),
                          #mask_sal = (['profile','level'], np.zeros((n_prof , n_obs_levels))*np.nan),
                          #mask_rho = (['profile','level'], np.zeros((n_prof , n_obs_levels))*np.nan),
                          obs_z =    (['profile','level'],    np.zeros((n_prof , n_obs_levels))*np.nan),
                          ex_mod_tem = (["profile", "ex_level"], mod_profiles.tem.values),
                          ex_mod_sal = (["profile", "ex_level"], mod_profiles.sal.values),
                          ex_depth =   (["profile", "ex_level"], mod_profiles.depth_0.values.T),
                          nn_ind_x = (["profile"], ind2D[0]),
                          nn_ind_y = (["profile"], ind2D[1])))

	# Vector of bools which save whether a datapoint is bad for some reason
    interp_dist = np.zeros(n_prof)
    
    # Now loop over profiles and interpolate model onto obs.
    for prof in range(0,n_prof):
        
        
        # Select the current profile
        mod_profile = mod_profiles.isel(profile = prof)
        obs_profile = en4.isel(profile = prof)
        
        # Check that model point is within threshold distance of obs
        # If not, skip vertical interpolation -> keep profile as nans
        interp_dist[prof] = coastgu.calculate_haversine_distance(
                                             obs_profile.longitude, 
                                             obs_profile.latitude, 
                                             mod_profile.longitude, 
                                             mod_profile.latitude)
        
        # Use bottom_level to mask dry depths
        if 'bottom_level' in mod_profile:
            bl = mod_profile.bottom_level.squeeze().values
            mod_profile = mod_profile.isel(z_dim=range(0,bl))
        
        # Interpolate model to obs depths using a linear interp
        # If interpolation fails for some or any reason, skip to next iteration        
        depth0 = mod_profile.depth_0.values
        depthnew = obs_profile.depth.values
                
        f = interp.interp1d(depth0, prof.mod_tem.values, fill_value=np.nan, 
                            bounds_error=False, kind = z_interp)
        data['mod_tem'][prof] = f( depthnew )
        f = interp.interp1d(depth0, prof.mod_sal.values, fill_value=np.nan, 
                            bounds_error=False, kind = z_interp)
        data['mod_sal'][prof] = f( depthnew )
        f = interp.interp1d(depth0, prof.obs_tem.values, fill_value=np.nan, 
                            bounds_error=False, kind = z_interp)
        data['obs_tem'][prof] = f( depthnew )
        f = interp.interp1d(depth0, prof.obs_sal.values, fill_value=np.nan, 
                            bounds_error=False, kind = z_interp)
        data['obs_sal'][prof] = f( depthnew )
        
        data['obs_z'][prof] = obs_profile.depth
        
        # # Calculate Density -- CURRENTLY NOT IMPLEMENTED ANYMORE
        # ap_obs = gsw.p_from_z( -obs_profile.depth, obs_profile.latitude )
        # ap_mod = gsw.p_from_z( -obs_profile.depth, mod_profile_int.latitude )
        # # Absolute Salinity            
        # sa_obs = gsw.SA_from_SP( obs_profile.sal, ap_obs, 
        #                         obs_profile.longitude, 
        #                         obs_profile.latitude )
        # sa_mod = gsw.SA_from_SP( mod_profile_int.sal, ap_mod, 
        #                         mod_profile_int.longitude, 
        #                         mod_profile_int.latitude )
        # # Conservative Temperature
        # ct_obs = gsw.CT_from_pt( sa_obs, obs_profile.tem ) 
        # ct_mod = gsw.CT_from_pt( sa_mod, mod_profile_int.tem ) 
        
        # # In-situ density
        # obs_rho = gsw.rho( sa_obs, ct_obs, ap_obs )
        # mod_rho = gsw.rho( sa_mod, ct_mod, ap_mod ) 
        
        # # Potential Density
        # obs_s0 = gsw.sigma0(sa_obs, ct_obs)
        # mod_s0 = gsw.sigma0(sa_mod, ct_mod)
        
        # Assign to main array
        #data['mod_rho'][prof] = mod_rho
        #data['obs_rho'][prof] = obs_rho
        #data['mod_s0'][prof] = mod_s0
        #data['obs_s0'][prof] = obs_s0
        
    print(' Interpolated Profiles.', flush=True)
    
    # Define seasons as month numbers and identify datapoint seasons
    month_season_dict = {1:0, 2:0, 3:2, 4:2, 5:2, 6:1,
                         7:1, 8:1, 9:3, 10:3, 11:3, 12:0}
    pd_time = pd.to_datetime(data.time.values)
    pd_month = pd_time.month
    season_save = [month_season_dict[ii] for ii in pd_month]
    
    data['season'] = ('profile', season_save)
    data.attrs['run_name'] = run_name
    data['interp_dist'] = ('profile', interp_dist)
    
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

def concatenate_output_files(files, fn_out):
    '''
    Uses xarray to concatenate a list of glob of files into a new single
    file.
    '''
    fn_list = glob.glob(files)
    ds_list = [xr.open_dataset(ff, chunks={'profile':10000}) for ff in fn_list]
    ds_concat = xr.concat(ds_list, dim='profile')
    write_ds_to_file(ds_concat, fn_out)
    return