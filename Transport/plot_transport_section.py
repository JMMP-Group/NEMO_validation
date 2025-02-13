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
    ds = xr.open_dataset("/gws/nopw/j04/jmmp/public/AMM15/DOMAIN_CFG/GEG_SF12.nc")
    ds = ds.isel(t=0)
    ds.to_netcdf("GEG_SF12.nc")
   

def _get_montly_means():
    """
    get monthly means of daily files
    """

    print ("n")
    years = np.arange(2004,2006,1)
    months = np.arange(1,13,1)
    start_date = "2004-01"
    end_date = "2006-01"
    dates = np.arange(start_date, end_date, dtype='datetime64[M]')
    t_series = []
    u_series = []
    v_series = []
    for date in dates:
        print (str(date))
        try:
            date_str = str(date).replace("-","")
            file_t=cfg.dn_dat + date_str + "01T0000Z_daily_grid_T.nc"
            file_u=cfg.dn_dat + date_str + "01T0000Z_daily_grid_U.nc"
            file_v=cfg.dn_dat + date_str + "01T0000Z_daily_grid_V.nc"
            chunks={"time_counter":1}
            t = xr.open_dataset(file_t, chunks=chunk).mean("time_counter")
            u = xr.open_dataset(file_u, chunks=chunk).mean("time_counter")
            v = xr.open_dataset(file_v, chunks=chunk).mean("time_counter")

            t = t.expand_dims(time_counter=[date])
            u = u.expand_dims(time_counter=[date])
            v = v.expand_dims(time_counter=[date])

            t_series.append(t)
            u_series.append(u)
            v_series.append(v)
        except Exception as e:
            print ("error: ", e)

    t_full = xr.concat(t_series, dim="time_counter")
    u_full = xr.concat(u_series, dim="time_counter")
    v_full = xr.concat(v_series, dim="time_counter")

    # save
    with ProgressBar():
        t_full.to_netcdf(cfg.dn_dat + "monthly_mean_grid_T.nc")
        u_full.to_netcdf(cfg.dn_dat + "monthly_mean_grid_U.nc")
        v_full.to_netcdf(cfg.dn_dat + "monthly_mean_grid_V.nc")

_get_montly_means()

def _get_flux():
    """
    use straitflux to get transport
    """

    help(master.transports)
    lon, lat = _get_ellet_line_positions()
    print (lon.data)
    model='CO9'
    product = 'volume'
    strait='Ellet' 
    time_start='2004-01'
    time_end='2014-12'
    path='Examples/' # path to save data

    #file_zv="GEG_SF12.nc"

    years = np.arange(2004,2014,1)
    for i in years:
        time_start=str(i)+'-01'
        time_end=str(i)+'-12'
        print(time_start,time_end)
        file_t=cfg.dn_dat + str(i) + "*01T0000Z_daily_grid_T.nc"
        file_u=cfg.dn_dat + str(i) + "*01T0000Z_daily_grid_U.nc"
        file_v=cfg.dn_dat + str(i) + "*01T0000Z_daily_grid_V.nc"
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
#_get_flux()

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
    
