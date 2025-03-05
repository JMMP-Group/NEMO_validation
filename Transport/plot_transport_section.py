import matplotlib as plt
from StraitFlux import masterscript_line as master
import xarray as xr
from PythonEnvCfg.config import config
cfg = config() # initialise variables in python
import numpy as np
from dask.diagnostics import ProgressBar


def _get_ellet_line_positions():
    """
    retrieve ellet line lat lon positions
    """

    path = cfg.dn_out + "transport/obs_for_ellet_line.nc"
    ds = xr.open_dataset(path)
    ds = ds.where(ds.time==2005, drop=True)

    return ds.longitude, ds.latitude

def remove_time_from_GEG_12():
    path = "/gws/nopw/j04/jmmp/public/AMM15/DOMAIN_CFG/GEG_SF12.nc"
    ds = xr.open_dataset(path)
    ds = ds.isel(t=0)
    ds.to_netcdf("GEG_SF12.nc")

def find_lat_lon_indicies(ds):
    """
    Find indicies using lat lon positions
    Saves compute resource
    """
   
    east = -8
    west = -14
    north = 58
    south = 56.5
    
    ind_w = min(np.abs(ds.nav_lon_grid_T - west).argmin("x_grid_T")).values
    ind_e = max(np.abs(ds.nav_lon_grid_T - east).argmin("x_grid_T")).values
    ind_n = max(np.abs(ds.nav_lat_grid_T - north).argmin("y_grid_T")).values
    ind_s = min(np.abs(ds.nav_lat_grid_T - south).argmin("y_grid_T")).values
    
    print ('n',ind_n)
    print ('s',ind_s)
    print ('e',ind_e)
    print ('w',ind_w)
   

def _get_montly_mean():
    """
    get monthly means of daily model files
    """

    def mean(start_date, end_date, vec="T"):

        if vec == "T":
            grid_var = f"_grid_{vec}"
        else:
            grid_var = ""
        dates = np.arange(start_date, end_date, dtype='datetime64[M]')
        ds_series = []
        for date in dates:
            try:
                date_str = str(date).replace("-","")
                fn=cfg.dn_dat + date_str + f"01T0000Z_daily_grid_{vec}.nc"
                chunks="auto"
                ds = xr.open_dataset(fn, chunks=chunks, decode_cf=True,
                                    decode_times=False)#.mean("time_counter")
                #ds = ds.drop("deptht_bounds")

                #find_lat_lon_indicies(ds)

                #dt = np.datetime64(date, "ns")
                #ds = ds.expand_dims(time_counter=[dt])

                n = 1064
                s = 837
                e = 620
                w = 197

                
                ds = ds.isel({f"x{grid_var}":slice(w,e),
                              f"y{grid_var}":slice(s,n)})
                # save
                #with ProgressBar():
                #    fn = f"{date_str}_Ellet_region_grid_{vec}.nc"
                #    save_path = cfg.dn_out + "transport/Ellet_cutout/" + fn
                #    ds.to_netcdf(save_path)


                ds_series.append(ds)

            except Exception as e:
                print ("error: ", e)

        full_series = xr.concat(ds_series, dim="time_counter")
        #full_series.time_counter.encoding["units"] = "seconds since 1900-01-01"
        #full_series.time_counter.encoding["dtype"] = "float64"
        #full_series.time_counter.attrs["dtype"] = "datetime64[ns]"

        # save
        with ProgressBar():
            date_range = (start_date + "_" + end_date).replace("-","")
            fn = f"{date_range}_Ellet_region_grid_{vec}.nc"
            save_path = cfg.dn_out + "transport/Ellet_cutout/" + fn
            full_series.to_netcdf(save_path)


    start_date = "2013-01"
    end_date = "2014-01"

    mean(start_date, end_date, vec="T")

#_get_montly_mean()

def _get_flux():
    """
    use straitflux to get transport
    """

    #help(master.transports)
    lon, lat = _get_ellet_line_positions()
    print (lon.data)
    model='CO9'
    product = 'volume'
    strait='Ellet' 
    time_start='2004-01'
    time_end='2014-12'
    path='Examples/' # path to save data

    #file_zv="GEG_SF12.nc"

    years = np.arange(2004,2005,1)
    for i in years:
        time_start=str(i)+'-01'
        time_end=str(i)+'-12'
        print(time_start,time_end)
        path = cfg.dn_out + "transport/Ellet_cutout/"
        file_t= path + str(i) + "*Ellet_region_grid_T.nc"
        file_u= path + str(i) + "*Ellet_region_grid_U.nc"
        file_v= path + str(i) + "*Ellet_region_grid_V.nc"
        transport = master.transports(product,
                                  strait,
                                  model,
                                  time_start,
                                  time_end,
                                  file_u,
                                  file_v,
                                  file_t,
                                  file_z=file_t,
                                  file_zu=file_u,
                                  file_zv=file_v,
                                  path_save=path,
                                  path_indices=path,
                                  path_mesh=path,
                                  set_latlon=True,
                                  lon_p=lon,
                                  lat_p=lat,
                                  Arakawa="Arakawa-C",
                                  saving=False)
        transport = transport.resample(time="1M").mean()
        print ("")
        print ("here")
        print (transport)
        print ("here")
        print (sdkjfh)
        print ("here")
_get_flux()

def plot_elet_obs_summary():
    """
    four panel plot of ladcp and ctd measurements
    """

    # intialise plots
    fig, axs = plt.subplots(2,3, figsize=(6.5,4))
    plt.subplots_adjust()

    # access data
    obs = xr.open_dataset("")
    
    # render time series of vels
    
