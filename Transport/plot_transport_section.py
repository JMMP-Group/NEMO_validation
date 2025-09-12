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
import os
import coast

class transport(object):

    def __init__(self, comp=True):
        self.comp = comp
        out_path_cross = cfg.dn_out + "transport/CrossSection/"
        out_path_cut = cfg.dn_out + "transport/Ellet_cutout/"
        os.makedirs(out_path_cross, exist_ok=True)
        os.makedirs(out_path_cut, exist_ok=True)
        if comp:
            out_path_cross = cfg.comp_case["proc_data"] + "transport/CrossSection/"
            out_path_cut = cfg.comp_case["proc_data"] + "transport/Ellet_cutout/"
            os.makedirs(out_path_cross, exist_ok=True)
            os.makedirs(out_path_cut, exist_ok=True)

    def _get_ellet_line_positions(self, sec=None):
        """
        retrieve ellet line lat lon positions
        """
    
        path = cfg.dn_out + "transport/obs_for_ellet_line.nc"
        ds = xr.open_dataset(path)
        ds = ds.where(ds.time==2005, drop=True)

        if sec == 'west':
            ds = ds.where(ds.longitude < -11)
        if sec == 'east':
            ds = ds.where(ds.longitude > -11)
    
        return ds.longitude, ds.latitude
    
    def remove_time_from_GEG_12(self):
        path = "/gws/nopw/j04/jmmp/public/AMM15/DOMAIN_CFG/GEG_SF12.nc"
        ds = xr.open_dataset(path)
        ds = ds.isel(t=0)
        ds.to_netcdf("GEG_SF12.nc")
    
    def find_lat_lon_indicies(self, ds):
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
       
    
    def _get_monthly_mean(self, path_in, path_out):
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
                print (date)
                #try:
                date_str = str(date).replace("-","")
                fn=path_in + date_str + f"*_25hourm_grid_{vec}.nc"
                #fn=path + date_str + f"01T0000Z_*_grid_{vec}.nc"
                print (fn)
                chunks="auto"
                #ds = xr.open_dataset(fn, chunks=chunks, decode_cf=True,
                #                    decode_times=False)#.mean("time_counter")
                #ds = xr.open_dataset(fn, chunks=chunks, decode_cf=True,
                #                    decode_times=False)#.mean("time_counter")
                nemo = coast.Gridded(fn, cfg.comp_case["grid"], multiple=True,
                                     config=cfg.fn_cfg_nemo)
                ds = nemo.dataset

                #ds = ds.drop("deptht_bounds")
    
                #find_lat_lon_indicies(ds)
    
                #dt = np.datetime64(date, "ns")
                #ds = ds.expand_dims(time_counter=[dt])
    
                n = 1064
                s = 837
                e = 620
                w = 197
    
                
                #ds = ds.isel({f"x{grid_var}":slice(w,e),
                #              f"y{grid_var}":slice(s,n)})
                ds = ds.isel({"x_dim":slice(w,e),
                              "y_dim":slice(s,n)})
                # save
                #with ProgressBar()om:
                #    fn = f"{date_str}_Ellet_region_grid_{vec}.nc"
                #    save_path = cfg.dn_out + "transport/Ellet_cutout/" + fn
                #    ds.to_netcdf(save_path)
    
    
                ds_series.append(ds)
    
                #except Exception as e:
                #    print ("error: ", e)
    
            #full_series = xr.concat(ds_series, dim="time_counter")
            full_series = xr.concat(ds_series, dim="t_dim")
            #full_series.time_counter.encoding["units"] = "seconds since 1900-01-01"
            #full_series.time_counter.encoding["dtype"] = "float64"
            #full_series.time_counter.attrs["dtype"] = "datetime64[ns]"
    
            # save
            with ProgressBar():
                date_range = (start_date + "_" + end_date).replace("-","")
                fn = f"{date_range}_Ellet_region_grid_{vec}.nc"
                save_path = path_out + "transport/Ellet_cutout/" + fn
                full_series.to_netcdf(save_path)
    
        #start_date = "2006-01"
        #end_date = "2007-01"
    
        for year in range(2006,2007):
            start_date = f"{year}-01"
            end_date = f"{year+1}-01"
            print (start_date)
            print (end_date)
            mean(start_date, end_date, vec="U")
            mean(start_date, end_date, vec="V")
    
    def _get_transport_all(self):
        """
        use straitflux to get transport for model and comparitor
        """
    
        self._get_transport(model=cfg.case, path=cfg.dn_out)
        if self.comp:
            self._get_transport(model=cfg.comp_case["case"], 
                            path=cfg.comp_case["proc_data"])
    
    def _get_transport(self, model, path, y0=2006, y1=2007, sec="all"):
    
        lon, lat = self._get_ellet_line_positions(sec=sec)
        product = 'volume'
        strait='Ellet' 
        save_path = path + "transport/CrossSection/"
        path = path + "transport/Ellet_cutout/"
    
        years = np.arange(y0,y1,1)
        for i in years:
            time_start=str(i)+'-01'
            time_end=str(i)+'-12'
            print(time_start,time_end)
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
                                      path_save=save_path,
                                      path_indices=save_path,
                                      path_mesh=save_path,
                                      set_latlon=True,
                                      lon_p=lon,
                                      lat_p=lat,
                                      Arakawa="Arakawa-C",
                                      saving=False)
    
            # get transport statistics
            transport = transport[model]
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
                        f"_Ellet_transport_stats_{sec}.nc"
                transport_stats.to_netcdf(path)
    
    
    def _get_cross_section(self):
    
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
    
    #_get_cross_section()

    def _get_transport_coast_format(self):
        """
        calculate transport wiht COAsT Methods
        """

        path_in=cfg.comp_case["raw_data"]
        lons, lats = self._get_ellet_line_positions()
        pts = list(zip(lats.values,lons.values))

        start_date = "2006-01"
        end_date = "2007-01"
        dates = np.arange(start_date, end_date, dtype='datetime64[M]')
        ds_series = []

        #n = 1064
        #s = 837
        #e = 620
        #w = 197

        nemo_f = coast.Gridded(fn_domain=cfg.comp_case["grid"],
                               config=cfg.fn_cfg_nemo_f)
        #nemo_f.dataset = nemo_f.dataset.isel({"x_dim":slice(w,e),
        #              "y_dim":slice(s,n)})
        for date in dates:
            print (date)
            date_str = str(date).replace("-","")
            chunks="auto"
            fn=path_in + date_str + f"*_25hourm_grid_U.nc"
            nemo_u = coast.Gridded(fn, cfg.comp_case["grid"], multiple=True,
                                   config=cfg.fn_cfg_nemo_u)
            fn=path_in + date_str + f"*_25hourm_grid_V.nc"
            nemo_v = coast.Gridded(fn, cfg.comp_case["grid"], multiple=True,
                                   config=cfg.fn_cfg_nemo_v)


            #nemo_u.dataset = nemo_u.dataset.isel({"x_dim":slice(w,e),
            #              "y_dim":slice(s,n)})
            #nemo_v.dataset = nemo_v.dataset.isel({"x_dim":slice(w,e),
            #              "y_dim":slice(s,n)})

            with ProgressBar():
                nemo_u.dataset = nemo_u.dataset.mean("t_dim").load()
                nemo_v.dataset = nemo_v.dataset.mean("t_dim").load()
            
            # save
            for i in range(len(pts) - 1):
                tran_f = coast.TransectF(nemo_f, pts[i], pts[i+1])
                tran_f.calc_flow_across_transect(nemo_u, nemo_v)
                vol = tran_f.data_cross_tran_flow.normal_transports.sum("r_dim")
                vol = vol.expand_dims("pts")
                # RDP note to self needs to be mean of pts i and pts i+1
                vol = vol.assign_coords(longitude=("pts",lons.data),
                                        latitude=("pts", lats.data))
                print (vol)
                print (sdfkj)

    def plot_ellet_transport(self, rolling=None):
        """
        plot time series of ellet transport
        """
    
        # initialise plots
        cm = 1/2.54  # centimeters in inches
        fig, ax = plt.subplots(1, figsize=(12*cm,6*cm))
        plt.subplots_adjust(bottom=0.2, top=0.95, right=0.95, left=0.15)
    
        # access data
        mod = xr.open_mfdataset(cfg.dn_out + 
                            "transport/ModelTransportStats/*transport*")
    
        mod_start = mod.time.min()
        mod_end = mod.time.max()
    
        mod = mod.resample(time="1MS").asfreq()/1e6
    
        ax.fill_between(mod.time, mod.quant.sel(quantile=0.25),
                                  mod.quant.sel(quantile=0.75))
        ax.plot(mod.time, mod["mean"], c='red')
    
    
        # observations
        path = cfg.dn_out + "transport/obs_for_ellet_line.nc"
        obs = xr.open_dataset(path)
        date = []
        for year, year_ds in obs.groupby("time"):
            date.append(datetime.datetime(year, int(year_ds.Month), 1))
        vol = obs.volume_transport / 1e6
        
        v_mean = vol.mean()
        v_lower = vol.mean() - vol.std()
        v_upper = vol.mean() + vol.std()
        pos = mod_end - np.timedelta64(12, "W")
        vl = ax.vlines(pos, v_lower, v_upper,
                             color='k', transform=ax.transData,
                             lw=2)
        vl.set(capstyle="round")
        ax.scatter(pos, vol.mean(), c='g')

        ax.set_xlim(mod_start,mod_end)

        ax.set_ylabel("Volume Transport (Sv)")
        ax.set_xlabel("Date")
    
        plt.show()
        plt.savefig("Figs/ellet_transport_timeseries.png")
    
    def plot_ellet_model_transport_cross_section(self):
        """
        plot cross section of velocities through Rockall Trough
        """
        def preprocess(ds):
            ds["x"] = ds.x.round(5)
            return ds
    
        path = cfg.dn_out + "transport/CrossSection/"
        mod = xr.open_mfdataset(path + "*cross*", preprocess=preprocess)
        lons = xr.open_mfdataset(path + "T_proj_points_CO9Ellet.nc").lon
        lats = xr.open_mfdataset(path + "T_proj_points_CO9Ellet.nc").lat
                
    
        # initialise plots
        fig, axs = plt.subplots(2)
    
        # need to plot mod and obs on same frame but for now plot separate 
        path = cfg.dn_out + "transport/obs_for_ellet_line.nc"
        obs = xr.open_dataset(path)
    
        vmin = -0.5 
        vmax = 0.5
    
        obs = obs.where(obs.time==2006, drop=True).squeeze()
        month = int(obs.Month.data)
    
        mod = mod.sel(time="2006-10").mean("time")
        mod = mod.sel(depth=slice(0,2500))
        print (mod)
        
        axs[0].pcolor(obs.Refdist, -obs.depth, obs.ladcp_velocity.T,
                      vmin=vmin, vmax=vmax, shading="auto", cmap=plt.cm.RdBu)
    
    
        #axs[1].pcolor(mod.x/1000, -mod.depth, mod.uv,
        #                    vmin=vmin, vmax=vmax, cmap=plt.cm.RdBu)
        plt.savefig("ellet_cross_200610.png")
    
    #plot_ellet_model_transport_cross_section()
    
    def get_climate_variables(self):
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
    
    def plot_elet_obs_summary(self):
        """
        four panel plot of ladcp and ctd measurements
        """
    
        # intialise plots
        fig, axs = plt.subplots(2,3, figsize=(6.5,4))
        plt.subplots_adjust()
    
        # access data
        obs = xr.open_dataset("")
        
        # render time series of vels

        
if __name__ == "__main__":
    trans = transport()
    #trans._get_monthly_mean(path_in=cfg.comp_case["raw_data"],
    #                        path_out=cfg.comp_case["proc_data"])
    #trans._get_monthly_mean(path_in=cfg.dn_dat,
    #                        path_out=cfg.dn_out)
    trans._get_transport_coast_format()
    #for year in [2005,2006,2007,2008,2009,2010,2011,2012,2013]:
    #    print ("year: ", year)
    #    trans._get_transport(model=cfg.case, 
    #               path=cfg.dn_out, y0=year,y1=year + 1, sec="east")
    #    trans._get_transport(model=cfg.case, 
    #               path=cfg.dn_out, y0=year,y1=year + 1, sec="west")
    #trans.plot_ellet_transport()
