import matplotlib.pyplot as plt
import numpy as np
print (np.__version__)
from StraitFlux import masterscript_line as master
from StraitFlux import masterscript_cross as master_cross
import StraitFlux
print (StraitFlux.__file__)
import xarray as xr
from PythonEnvCfg.config import config
cfg = config() # initialise variables in python
from dask.diagnostics import ProgressBar
import datetime


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


    start_date = "2012-01"
    end_date = "2013-01"

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

        # get transport statistics
        transport = transport.CO9
        mean = transport.resample(time="1MS").mean()
        mean.name = "mean"
        quant = transport.resample(time="1MS").quantile([0.25,0.5,0.75])
        quant.name = "quant"
        std = transport.resample(time="1MS").std()
        std.name = "std"

        # merge into ds
        transport_stats = xr.merge([mean,quant,std])

        # save
        with ProgressBar():
            path = cfg.dn_out + "transport/" + str(i) + \
                    "_Ellet_transport_stats.nc"
            transport_stats.to_netcdf(path)
#_get_flux()

def _get_cross_section():

    lon, lat = _get_ellet_line_positions()
    print (lon)
    model='CO9'
    product = 'volume'
    strait='Ellet' 

    out_path = cfg.dn_out + "transport/CrossSection/"
    years = np.arange(2004,2005,1)
    for i in years:
        time_start=str(i)+'-01'
        time_end=str(i)+'-12'
        print(time_start,time_end)
        path = cfg.dn_out + "transport/Ellet_cutout/"
        file_t= path + str(i) + "*Ellet_region_grid_T.nc"
        file_u= path + str(i) + "*Ellet_region_grid_U.nc"
        file_v= path + str(i) + "*Ellet_region_grid_V.nc"

        uv=master_cross.vel_projection(strait,
                                  model,
                                  time_start,
                                  time_end,
                                  file_u,
                                  file_v,
                                  file_t,
                                  file_z=file_t,
                                  file_zu=file_u,
                                  file_zv=file_v,
                                  set_latlon=True,
                                  lon_p=lon,
                                  lat_p=lat,
                                  Arakawa="Arakawa-C",
                                  saving=True,
                                  path_save=out_path)


        # save
        #with ProgressBar():
        #    path = cfg.dn_out + "transport/CrossSection/" + str(i) + \
        #            "_Ellet_velocity_cross_section.nc"
        #    uv.to_netcdf(path)#, encoding={"time": {"dtype": "i4"}})

_get_cross_section()

def plot_ellet_transport():
    """
    plot time series of ellet transport
    """

    # initialise plots
    fig, ax1 = plt.subplots(1)

    ax2 = ax1.twinx()

    # access data
    print (cfg.dn_out + "transport/ModelTransportStats/*")
    mod = xr.open_mfdataset(cfg.dn_out + "transport/ModelTransportStats/*transport*")

    mod_start = mod.time.min()
    mod_end = mod.time.max()
    NAO = get_climate_variables()
    NAO = NAO.sel(time=slice("2006-02-01",mod_end))
    mod = mod.sel(time=slice("2006-02-01",mod_end))
    NAO = NAO.sel(time=slice(mod_start,mod_end))

    mod = mod.resample(time="1MS").asfreq()/1e6

    NAO = NAO.rolling(time=12).mean()
    mod = mod.rolling(time=12).mean()
    #mod = mod/abs(mod).max("time")
    #NAO = NAO/abs(NAO).max("time")
    print (NAO.max())
    print (mod.max())
    print (NAO.min())
    print (mod.min("time"))
    cov = np.ma.correlate(np.ma.masked_invalid(NAO), np.ma.masked_invalid(mod["mean"]), mode="same")
    print (NAO)
    print (mod["mean"])
    print (cov)
    #print (skjdfh)
    ax1.fill_between(mod.time, mod.quant.sel(quantile=0.25),
                              mod.quant.sel(quantile=0.75))
    ax1.plot(mod.time, mod["mean"], c='red')

    ax2.plot(NAO.time, NAO, c="orange")

    # observations
    path = cfg.dn_out + "transport/obs_for_ellet_line.nc"
    obs = xr.open_dataset(path)
    date = []
    for year, year_ds in obs.groupby("time"):
        print (year)
        print (year_ds)
        date.append(datetime.datetime(year, int(year_ds.Month), 1))
    vol = obs.volume_transport / 1e4
    print (date)
    plt.scatter(date, vol, c='g')

    plt.show()

def plot_ellet_model_transport_cross_section():
    """
    plot cross section of velocities through Rockall Trough
    """
    def preprocess(ds):
        ds["x"] = ds.x.round(5)
        return ds

    path = cfg.dn_out + "transport/CrossSection/*cross*"
    mod = xr.open_mfdataset(path, preprocess=preprocess)

    print (mod)

    # initialise plots
    fig, ax = plt.subplots(1)

    print (sjdk)



#plot_ellet_model_transport_cross_section()

def get_climate_variables():
    """
    get NAO
    """

    NAO = np.loadtxt(cfg.dn_out + "transport/NAO/nao_station_monthly.txt",
                     skiprows=2)
    NAO_data = NAO[:,1:].flatten()
    NAO_years = NAO[:,0].astype("int")
    NAO_time = np.arange(str(NAO_years[0]) + "-01",
                         str(NAO_years[-1]+1) + "-01",
                             dtype="datetime64[M]")
    print (NAO_time)
    

    NAO_xr = xr.DataArray(NAO_data, dims=("time"),
                          coords={"time": NAO_time})

    return NAO_xr


#plot_ellet_transport()

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
    
