from PythonEnvCfg.config import config, bounds
cfg = config() # initialise variables in python

import coast
import xarray as xr
import os
import copernicusmarine
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.interpolate import griddata
import numpy as np
import xeofs as xe
import glob
import time
from dask.diagnostics import ProgressBar
import matplotlib
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
import cartopy.feature as cfeature

matplotlib.rcParams.update({'font.size': 8})

class extract_surface(object):
    def __init__(self):

        # paths
        self.fn_dom = cfg.dn_dom + cfg.grid_nc
        self.fn_dat = cfg.dn_out + "profiles/gridded*.nc"
        self.fn_out = cfg.dn_out + 'surface_maps/'

    def surface_state_climatology_native_model(self):
        """ 
        Create surface climatology for full model data

        Open gridded model data and calculate climatology. The gridded model
        dataset is on the native model grid.
        """

        ds = xr.open_mfdataset(self.fn_dat, combine="nested", 
                               concat_dim="t_dim", parallel=True)
        clim = coast.Climatology()
        clim_mean = clim.make_climatology(ds, "season",
                       fn_out=self.fn_out + "surface_state_climatology.nc")

    def surface_state_climatology_binned_model(self):
        """
        Create surface (5 m) climatoligy from binned surface model data.
        """

        ds = xr.open_mfdataset(self.fn_out, combine="nested", 
                               concat_dim="t_dim", parallel=True)
        ds = ds[["temperature","salinity"]]
        clim = coast.Climatology()
        clim_mean = clim.make_climatology(ds, "season",
              fn_out=self.fn_out + "surface_state_climatology_binned_model.nc")

class satellite(object):
    """
    Class for validating against satellite data
    """

    def __init__(self, cfg_fn):
        self.cfg_fn = cfg_fn

    def create_cmems_login(self):
        """ create login config file """

        copernicusmarine.login()

    def get_cmems(self, var="ssh"):
        """ download cmems data """

        bdy = bounds("AMM15")

        data_request = {
           "longitude" : [bdy.lonbounds[0], bdy.lonbounds[1]],
           "latitude" : [bdy.latbounds[0], bdy.latbounds[1]],
           "time" : ["2004-01-01", "2014-01-01"],
        }
        
        if var == "sst":
           data_request["fn"] = "METOFFICE-GLO-SST-L4-REP-OBS-SST"
           self.var_str = "analysed_sst"

        if var == "ssh":
           data_request["fn"] = "cmems_obs-sl_eur_phy-ssh_my_allsat-l4-duacs-0.0625deg_P1D"
           self.var_str = "adt"
        if var == "uvel":
           data_request["fn"] = "cmems_obs-sl_eur_phy-ssh_my_allsat-l4-duacs-0.0625deg_P1D"
           self.var_str = "ugos"
        if var == "vvel":
           data_request["fn"] = "cmems_obs-sl_eur_phy-ssh_my_allsat-l4-duacs-0.0625deg_P1D"
           self.var_str = "vgos"

        data_request["variables"] = [self.var_str]
        # Load xarray dataset
        self.ds = copernicusmarine.open_dataset(
            dataset_id = data_request["fn"],
            minimum_longitude = data_request["longitude"][0],
            maximum_longitude = data_request["longitude"][1],
            minimum_latitude = data_request["latitude"][0],
            maximum_latitude = data_request["latitude"][1],
            start_datetime = data_request["time"][0],
            end_datetime = data_request["time"][1],
            variables = data_request["variables"]
        )[self.var_str]

        return self.ds

    def get_l4_SWOT(self):
        """ get geostrophic currents from compilation that incl. SWOT """

        self.ds = xr.open_mfdataset(cfg.dn_swot + "*.nc", chunks="auto")

    def monthly_mean(self):
        """ average over month """

        self.ds = self.ds.resample(time='1MS').mean().load()

    def interpolate_to_model(self, src_da):
        """ interpolate lat-lon to horizontal grid """

        domcfg = xr.open_dataset(self.cfg_fn)

        tgt_lon =  domcfg.nav_lon
        tgt_lat =  domcfg.nav_lat
        
        target = (tgt_lon, tgt_lat)

        src_lon = src_da.longitude.values
        src_lat = src_da.latitude.values

        src_mlon, src_mlat = np.meshgrid(src_lon, src_lat)

        points = (src_mlon.flatten(), src_mlat.flatten())

        n_grid_all = []
        for i, (time, da_t) in enumerate(src_da.groupby("time")):
            print (time)
            values = da_t.values.flatten()
            values_masked = (np.nan_to_num(values, nan=-9999))
            
            if i == 0:
                n_grid_all = griddata(points, values_masked, target,
                                      method="nearest")[:,:,np.newaxis]
            else:
                n_grid = griddata(points, values_masked, target,
                          method="nearest")[:,:,np.newaxis]
                n_grid_all = np.concatenate([n_grid_all, n_grid], axis=2)

        n_grid_all = xr.where(n_grid_all == -9999, np.nan, n_grid_all)

        src_da = xr.DataArray(
                             data=n_grid_all,
                             dims=["y","x","time"],
                             coords={"longitude": (["y","x"],tgt_lon.values),
                                     "latitude": (["y","x"],tgt_lat.values),
                                     "time": src_da.time},
                             name=self.var_str)#.to_dataset()
        return src_da

    def quick_compare(self):
        """ quick plot """

        snapshot_obs = self.ds.isel(time=0).adt
        snapshot_co9 = self.ds.isel(time=0)

        plt.pcolor(snapshot)
        plt.show()

    def save_ds(self, fn_name):
        """ save satellite data """

        fn = f"{cfg.dn_out}/satellite/{fn_name}"
        self.ds.to_netcdf(fn)

    def get_geostrophic_vels_CMEMS(self, time_step):
        """ get CMEMS data """

        # get surface geostrophic velocities
        path = f"{cfg.dn_out}/satellite/"
        u_fn = path + f"CMEMS_L4_satellite_uvel_{time_step}_gridded_to_P2.0.nc"
        v_fn = path + f"CMEMS_L4_satellite_vvel_{time_step}_gridded_to_P2.0.nc"
        u = xr.open_dataarray(u_fn, chunks="auto")
        v = xr.open_dataarray(v_fn, chunks="auto")

        return u, v

    def get_geostrophic_vels_AVISO_SWOT(self, time_step):
        """ get AVISO SWOT data """

        path = f"{cfg.dn_out}/satellite/"
        fn = path + f"dt_allsat_phy_l4_{time_step}_gridded_to_P2.0.nc"
        aviso_swot = xr.open_dataset(fn)
        u = aviso_swot.ugos
        v = aviso_swot.vgos

        return u, v


    def get_KE(self, save=False, time_step="1MS", src="CMEMS"):
        """ get mean and eddy kinetic energy of surface currents """

        if src == "CMEMS":
            u, v = self.get_geostrophic_vels_CMEMS(time_step)
        elif src == "AVISO_SWOT":
            u, v = self.get_geostrophic_vels_AVISO_SWOT(time_step)

        # time mean
        u_bar = u.mean("time")
        v_bar = v.mean("time")

        # devaition from mean
        u_prime = u_bar - u
        v_prime = v_bar - v

        # calculate mean and eddy kinetic energy
        MKE = 0.5 * (u_bar**2 + v_bar**2)
        EKE = 0.5 * ((u_prime**2).mean("time") +
                     (v_prime**2).mean("time"))

        # label vars
        MKE.name = "MKE"
        EKE.name = "EKE"

        # merge into single dataset
        KE = xr.merge([MKE,EKE])

        # save
        if save: 
            d0 = u.time.min().dt.date.values
            d1 = u.time.max().dt.date.values

            fn_str = f"{d0}_{d1}_{time_step}_{src}_satellite_KE.nc"
            fn = f"{cfg.dn_out}/satellite/" + fn_str
            KE.to_netcdf(fn)

        return KE

class model_surface(object):

    def __init__(self, fn_raw, fn_proc, src_t_coord="time", src_x_coord="x",
                                        src_y_coord="y", src_z_coord="z"):
        self.fn_path = fn_raw
        self.fn_proc = fn_proc
        self.src_t_coord = src_t_coord
        self.src_x_coord = src_x_coord
        self.src_y_coord = src_y_coord
        self.src_z_coord = src_z_coord

        eof_norm=False
        eof_std=True

    def map_dimension_coords(self, ds):
 
        """
        create uniform coordinate variables

        note: might better be handled by coast

        """
        
        print (ds)
        ds = ds.rename({self.src_t_coord:"time",
                        self.src_x_coord:"x",
                        self.src_y_coord:"y",
                        self.src_z_coord:"z"})

        return ds

    def map_lat_lon_names(self, src_nav_lon, src_nav_lat):
        """
        rename nav_lat and nav_lon
        """

        self.ds = self.ds.rename({src_nav_lon:"longitude",
                                  src_nav_lat:"latitude"})

    def get_domain_cfg(self, rename=None):

        domcfg = xr.load_dataset(self.cfg_fn).squeeze()

        if rename:
            domcfg = domcfg.rename({rename:"bathy"})

        return domcfg

    def get_time_mean_var(self, freq="1MS", grid="T", var_nam="sossheig",
                          save=True):
        drange = np.arange(f"{cfg.y0}-01", f"{cfg.y1}-01",
                    dtype="datetime64[M]")
        yrange = np.arange(int(cfg.y0), int(cfg.y1))
        t0 = time.time()
        chunks = {"time_counter":1}
        path_list = []
        for y in yrange:
            paths = glob.glob(self.fn_path + f"{y}*_25hourm_grid_{grid}.nc")
            path_list += paths

        da_var = xr.open_dataset(path_list[0], chunks=chunks)[var_nam]
        for path in path_list[1:]:
            da = xr.open_dataset(path, chunks=chunks)[var_nam]
            da_var = xr.concat([da_var, da], dim="time_counter")

        da_var = self.map_dimension_coords(da_var)

        # check for depth var, flawed if not spatial/time dims present
        if len(da_var.dims) > 3:
            da_var = da_var.isel(z=0)

        if freq:
            da_var = da_var.resample(time=freq).mean()

        with ProgressBar():
            self.ds = da_var.load()
        t1 = time.time()
        print ((t1-t0)/60)

        if save:
            if freq == None:
                freq="25hourm"
            self.ds.to_netcdf(self.fn_proc +
                             f"{cfg.y0}_{cfg.y1}_{freq}_{var_nam}.nc")

    def interpolate_sp_to_model(self, tgt, src):
        """ interpolate lat-lon to horizontal grid """

        tgt_lon = tgt.longitude
        tgt_lat = tgt.latitude
        
        target = (tgt_lon, tgt_lat)

        src_lon = src.longitude
        src_lat = src.latitude

        # check longitude format
        if src_lon.max() > 180:
            src_lon = xr.where(src_lon > 180, src_lon - 360, src_lon)

        if len(src_lon.shape) == 1:
            src_lon, src_lat = np.meshgrid(src_lon, src_lat)
            points = (src_lon.flatten(), src_lat.flatten())
        else:

            points = (src_lon.data.flatten(), src_lat.data.flatten())

        n_grid = []
        for time, ds_t in src.groupby("time"):
            values = (ds_t.data.flatten())
            print ("values", values.shape)
            
            n_grid.append(
             griddata(points, values, target, method="cubic")[:,:,np.newaxis])

        n_grid_all = np.concatenate(n_grid, axis=2)

        ds = xr.DataArray(
                             data=n_grid_all,
                             dims=["y","x","time"],
                             coords={"longitude": (["y","x"],tgt_lon.values),
                                     "latitude": (["y","x"],tgt_lat.values),
                                     "time": src.time},
                             name="sp")
        return ds

    def remove_inverse_barometer(self, src="era5"):
        """
        UNDER CONSTRUCTION
        remove atmospheric loading 
        """

        g = 9.80665
        rho = 1026
        g_rho = g * rho
        pref = 101000

        # getting this variable is tricky
        # it requires the raw surface foring being interpolated
        # with the same interpolation method use during model run
        # apr = ...

        #ssh_ib = - ( apr - pref ) / g_rho 

        if src == "era5":
            print (cfg.dn_era5)
            fn_list = [f"{cfg.dn_era5}ERA5_sp_y{y}.nc" for y in
                       range(int(cfg.y0),int(cfg.y1))]
            ds_list = []
            for fn in fn_list:
                print (fn)
                sp_year = xr.open_dataarray(fn, chunks="auto")
                sp_mean = sp_year.resample(time="MS").mean("time")
                ds_list.append(sp_mean.load())
            sp = xr.concat(ds_list, "time")
        if src == "era_interim":

            def preprocess(ds):
                ds = ds.expand_dims(dict(time=[ds.time.data]))
                return ds
            path_list = []
            for y in range(int(cfg.y0),int(cfg.y1)):
                paths = glob.glob(f"{cfg.dn_era_interim}ei.moda.an.sfc.regn128sc.{y}*.grib")
                path_list += paths
            #path = cfg.dn_era_interim+ "ei.moda.an.sfc.regn128sc.*100.grib"
            #           for y in range(int(cfg.y0),int(cfg.y1))]
            sp_list = []
            for fn in path_list:
                print (fn)
                sp_year = xr.open_dataset(fn, engine="cfgrib"
                           ).sp
                sp_year = sp_year.expand_dims(dict(time=[sp_year.time.data]))
                sp_list.append(sp_year)
            sp = xr.concat(sp_list, "time")

        # interpolate to model grid
        sp = self.interpolate_sp_to_model(self.ds, sp)

        ssh_ib = (sp - pref) / g_rho

        self.ds = self.ds + ssh_ib

    def get_KE(self, u, v):
        """ get mean and eddy kinetic energy of surface currents """

        uT = 0.5 * (u + u.shift(y=1))
        vT = 0.5 * (v + v.shift(y=1))

        with ProgressBar():
            uT_bar = uT.mean(self.src_t_coord).load()
            vT_bar = vT.mean(self.src_t_coord).load()

        uT_prime = uT_bar - uT
        vT_prime = vT_bar - vT

        MKE = 0.5 * (uT_bar**2 + vT_bar**2)
        with ProgressBar():
            EKE = 0.5 * ((uT_prime**2).mean(self.src_t_coord) +
                         (vT_prime**2).mean(self.src_t_coord)).load()

        MKE.name = "MKE"
        EKE.name = "EKE"

        KE = xr.merge([MKE,EKE])

        return KE

class satellite_plot(object):

    def plot_model_and_satellite_snapshot_ssh(self, mod_ssh, sat_ssh):
        """ plot ssh for model and satellite """

        fig, axs = plt.subplots(2)

        # time slice
        mod_ssh = mod_ssh.sel(time_counter="2004-01")
        sat_ssh = sat_ssh.sel(time="2004-01")

        p0 = axs[0].pcolor(mod_ssh.squeeze(), vmin=-1, vmax=1)
        p1 = axs[1].pcolor(sat_ssh.squeeze(), vmin=-1, vmax=1)

        plt.colorbar(p0, ax=axs[0])
        plt.colorbar(p1, ax=axs[1])
        plt.show()

    def get_data_lims(self, da_list):
        """ get absolute maximim of data """

        lims = []
        for  da in da_list:
            lims.append(abs(da).max().values)

        lim = max(lims)

        return -lim, lim

    def plot_eof_validation(self, mod0, mod1, sat, var):
        """ plot eof breakdown of model versus obs """

        # initialise figure
        fig = plt.figure(figsize=(6.5,6.5))

        # initialise gridspec
        gs0 = gridspec.GridSpec(ncols=2, nrows=1)
        gs1 = gridspec.GridSpec(ncols=2, nrows=3)

        ## set frame bounds
        gs0.update(top=0.96, bottom=0.84, left=0.12, wspace=0.6, hspace=0.12,
                   right=0.85)
        gs1.update(top=0.77, bottom=0.07, left=0.12, wspace=0.6, hspace=0.08,
                   right=0.85)

        # get domain cfg
        domcfg = xr.open_dataset(cfg.dn_dom + cfg.grid_nc,
                                 chunks="auto").squeeze()

        # set projection
        mid_lat = np.mean([44,domcfg.nav_lat.max().values])
        mid_lon = -2.6
        proj=ccrs.EquidistantConic(central_latitude=mid_lat,
          standard_parallels=(44,domcfg.nav_lat.max().values),
          central_longitude=mid_lon) 
        plt_proj=ccrs.PlateCarree()

        # assign axes to lists
        axs0 = []
        for i in range(2):
            axs0.append(fig.add_subplot(gs0[i]))
        axs1 = []
        for i in range(6):
            axs1.append(fig.add_subplot(gs1[i], projection=proj))

        # get data
        path = f"{cfg.dn_out}/satellite/"
        eof_type = ""
        mod0_scores = xr.open_dataarray(
                           f"{path}{mod0}_{var}_eof{eof_type}_scores.nc")
        mod0_comp = xr.open_dataarray(
                           f"{path}{mod0}_{var}_eof{eof_type}_components.nc")

        sat_scores = xr.open_dataarray(
                           f"{path}{sat}_{var}_eof{eof_type}_scores.nc")
        sat_comp = xr.open_dataarray(
                           f"{path}{sat}_{var}_eof{eof_type}_components.nc")

        #path = cfg.comp_case["proc_data"] + "/satellite/"
        mod1_scores = xr.open_dataarray(
                           f"{path}{mod0}_{var}_eof{eof_type}_scores.nc")
        mod1_comp = xr.open_dataarray(
                           f"{path}{mod0}_{var}_eof{eof_type}_components.nc")

        # get p-value
        mod0_pval = xr.corr(mod0_scores, sat_scores, dim="time")
        mod1_pval = xr.corr(mod1_scores, sat_scores, dim="time")

        # get score lims
        scores_list_m1 = [mod0_scores.sel(mode=1),
                          mod1_scores.sel(mode=1),
                          sat_scores.sel(mode=1)]
        scores_list_m2 = [mod0_scores.sel(mode=2),
                          mod1_scores.sel(mode=2),
                          sat_scores.sel(mode=2)]

        self.vmin_m1_s, self.vmax_m1_s =  self.get_data_lims(scores_list_m1)
        self.vmin_m2_s, self.vmax_m2_s =  self.get_data_lims(scores_list_m2)

        # get component lims
        comp_list_m1 = [mod0_comp.sel(mode=1),
                        mod1_comp.sel(mode=1),
                        sat_comp.sel(mode=1)]
        comp_list_m2 = [mod0_comp.sel(mode=2),
                        mod1_comp.sel(mode=2),
                        sat_comp.sel(mode=2)]

        self.vmin_m1_c, self.vmax_m1_c =  self.get_data_lims(comp_list_m1)
        self.vmin_m2_c, self.vmax_m2_c =  self.get_data_lims(comp_list_m2)

        def render(axs, comp, scores, i, label, plt_proj):
            axs0[0].plot(scores.time, scores.sel(mode=1), label=label,
                          lw=1.0)
            axs0[1].plot(scores.time, scores.sel(mode=2), label=label,
                          lw=1.0)

            #for ax in axs0:
            #    print (scores.time)
            #    n_steps = np.ceil(scores.time.max() - scores.time.min() / 5)
            #    print (n_steps)
            #    print (nsdflkj)
            #    ax.set_xticks(np.arange)

            
            p0 = axs1[i*2].pcolormesh(comp.longitude, comp.latitude,
                               comp.sel(mode=1).squeeze(), transform=plt_proj,
                               vmin=self.vmin_m1_c, vmax=self.vmax_m1_c,
                               cmap=plt.cm.RdBu_r)
            p1 = axs1[i*2+1].pcolormesh(comp.longitude, comp.latitude,
                               comp.sel(mode=2).squeeze(), transform=plt_proj,
                               vmin=self.vmin_m2_c, vmax=self.vmax_m2_c,
                               cmap=plt.cm.RdBu_r)
            #for ax in axs1[i*2:i*2+2]:
            #    ax.pcolormesh(comp.longitude, comp.latitude,
             #        xr.where(comp.isel(mode=1).squeeze() == np.nan, 1, np.nan))


            # set extent
            #lon0 = comp.longitude.isel(x=0, y=0).values
            lon0 = -15
            lon1 = 9.8
            axs[i*2].set_extent([lon0, lon1, 46, 62])
            axs[i*2+1].set_extent([lon0, lon1, 46, 62])
            #axs1[i].set_xlim(lon0, lon1)
            #axs1[i].set_ylim(46, 62)
            #axs1[i+3].set_xlim(lon0, lon1)
            #axs1[i+3].set_ylim(46, 62)
            return p0, p1

        render(axs1, mod0_comp, mod0_scores, 0, "CO9", plt_proj)
        render(axs1, mod1_comp, mod1_scores, 1, "CO7", plt_proj)
        print (sat_comp)
        #sat_comp = sat_comp.rename({"nav_lon":"longitude",
        #                            "nav_lat":"latitude"})
        p0, p1 = render(axs1, sat_comp, sat_scores, 2, "Sat. L4", plt_proj)

        axs0[1].legend(loc="upper left", bbox_to_anchor=(1.02,1),
                        bbox_transform=axs0[1].transAxes)

        # set timeseries lims and ticks
        for ax in axs0:
            ticks = ax.get_xticks()
            tick_labels = ax.get_xticklabels()
            ax.set_xlim(mod0_scores.time.min(), mod0_scores.time.max())
            
            year_min = mod0_scores.time.min().dt.year.values
            year_max = mod0_scores.time.max().dt.year.values
            step = np.ceil(
                     (year_max - year_min)/ 5)
            ticks = np.arange(str(year_min+1), str(year_max+1), int(step),
                                dtype="datetime64[Y]")
            ax.set_xticks(ticks)
            ax.set_xticklabels(ticks)

            # add zero line
            ax.axhline(0, ls="--", lw=0.5, c="k", zorder=10)
            

        axs0[0].set_ylim(self.vmin_m1_s*1.05,self.vmax_m1_s*1.05)
        axs0[1].set_ylim(self.vmin_m2_s*1.05,self.vmax_m2_s*1.05)

        # add p-vals
        #p =  str(np.round(mod0_pval.sel(mode=1).data, 2))
        #axs0[0].text(0.5, 0.95, "p = " + p,
        #              ha="left", va="top", transform=axs0[0].transAxes)
        #p =  str(np.round(mod1_pval.sel(mode=1).data, 2))
        #axs0[0].text(0.75, 0.95, "p = " + p,
        #              ha="left", va="top", transform=axs0[0].transAxes)
        #p =  str(np.round(mod0_pval.sel(mode=2).data, 2))
        #axs0[1].text(0.5, 0.95, "p = " + p,
        #              ha="left", va="top", transform=axs0[1].transAxes)
        #p =  str(np.round(mod1_pval.sel(mode=2).data, 2))
        #axs0[1].text(0.75, 0.95, "p = " + p,
        #              ha="left", va="top", transform=axs0[1].transAxes)

        # colorbar
        pos0 = axs1[0].get_position()
        pos1 = axs1[2].get_position()
        pos2 = axs1[4].get_position()

        cbar_ax = fig.add_axes([pos0.x1+0.02, (pos2.y1+pos2.y0)/2, 
                              0.02, (pos0.y1+pos0.y0)/2 - (pos2.y1+pos2.y0)/2])

        cbar = fig.colorbar(p0, cax=cbar_ax, orientation='vertical')
        cbar.ax.text(4.2, 0.5, f"{var.upper()} EOF Compontent",
                     rotation=90, transform=cbar.ax.transAxes,
                     va='center', ha='left')

        pos0 = axs1[1].get_position()
        pos1 = axs1[3].get_position()
        pos2 = axs1[5].get_position()

        #cbar_ax = fig.add_axes([0.85, 0.18, 
        #                        0.02, 0.47])
        cbar_ax = fig.add_axes([pos0.x1+0.02, (pos2.y1+pos2.y0)/2, 
                              0.02, (pos0.y1+pos0.y0)/2 - (pos2.y1+pos2.y0)/2])

        cbar = fig.colorbar(p1, cax=cbar_ax, orientation='vertical')
        cbar.ax.text(4.2, 0.5, f"{var.upper()} EOF Component",
                     rotation=90, transform=cbar.ax.transAxes,
                     va='center', ha='left')

        # set labels
        axs0[0].set_title("Mode 1")
        axs0[1].set_title("Mode 2")

        axs1[0].text(0.05,0.95, "CO9", ha="left", va="top",
                      transform=axs1[0].transAxes)
        axs1[2].text(0.05,0.95, "CO7", ha="left", va="top",
                      transform=axs1[2].transAxes)
        axs1[4].text(0.05,0.95, "Sat. L4", ha="left", va="top",
                      transform=axs1[4].transAxes)

        axs1[1].text(0.05,0.95, "CO9", ha="left", va="top",
                      transform=axs1[1].transAxes)
        axs1[3].text(0.05,0.95, "CO7", ha="left", va="top",
                      transform=axs1[3].transAxes)
        axs1[5].text(0.05,0.95, "Sat. L4", ha="left", va="top",
                      transform=axs1[5].transAxes)

        for ax in axs1[:4]:
            ax.set_xticklabels([])
        for ax in axs1[1::2]:
            ax.set_yticklabels([])
        for ax in axs1[::2]:
            ax.set_ylabel("Latitude")
        for ax in axs1[4:]:
            ax.set_xlabel("Longitude")
        for ax in axs0:
            ax.set_xlabel("Year")
            ax.set_ylabel("EOF Score")

        #cfg_fn = '/gws/nopw/j04/jmmp/public/AMM15/DOMAIN_CFG/GEG_SF12.nc'
       # domcfg = xr.open_dataset(cfg_fn).squeeze()
       # print (domcfg)
       # msk = xr.where(domcfg.top_level == 0, 0, np.nan)
       # print (msk)
        for i, ax in enumerate(axs1):
            if i in [4,5]: 
                b_label=True
            else:
                b_label=False

            self.add_land_and_gridlines(ax, domcfg, plt_proj, r_label=False,
                           t_label=False, b_label=b_label, l_label=True)
            ##ax.pcolormesh(domcfg.nav_lon, domcfg.nav_lat, msk, cmap="grey")
            ##ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')
            #land_50m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
            #                            edgecolor='k',
            #                            facecolor='grey')
            #ax.add_feature(land_50m, zorder=100, lw=0.5)

            #ax.set_xticks([-15, -10, -5, 0, 5, 10], crs=ccrs.PlateCarree())
            #ax.set_yticks([50, 55, 60], crs=ccrs.PlateCarree())
            #lon_formatter = LongitudeFormatter(zero_direction_label=True)
            #lat_formatter = LatitudeFormatter()

        # set save path
        eof_str = eof_type.lstrip("_")
        if eof_str == "":
            fig_path = f"FIGS/CO9_CO7_CMEMS_Satellite_{var}_pca.png"
        else:
            fig_path = f"FIGS/CO9_CO7_CMEMS_Satellite_{var}_{eof_str}_pca.png"
        print (fig_path)

        # save figure
        plt.savefig(fig_path, dpi=600)

    def plot_KE_by_component(self, KE_var="MKE", vmin=0, vmax=0.06):
        """
        plot KE for satellite data and two models
        a 2x3 panel plot with:
           - row 1 : monthly EKE for all cases
           - row 2 : daily EKE for models
        """
        # initialise figure
        fig, axs = plt.subplots(2, 3, figsize=(5.5,6.5))

        # get primary model data
        path = f"{cfg.dn_out}/satellite/"
        p_mod_d = xr.open_dataset(path + "2004_2014_25hourm_KE.nc",
                                  chunks="auto")[KE_var]
        p_mod_m = xr.open_dataset(path + "2004_2014_monthly_KE.nc",
                                  chunks="auto")[KE_var]

        # get satellite data
        sat_d = xr.open_dataset(path + "2004_2014_1D_satellite_KE.nc",
                                  chunks="auto")[KE_var]
        sat_m = xr.open_dataset(path + "2004_2014_1MS_satellite_KE.nc",
                                  chunks="auto")[KE_var]

        # get comparison model data
        path = cfg.comp_case["proc_data"] + "/satellite/"
        c_mod_d = xr.open_dataset(path + "2004_2014_25hourm_KE.nc",
                                  chunks="auto")[KE_var]
        c_mod_m = xr.open_dataset(path + "2004_2014_monthly_KE.nc",
                                  chunks="auto")[KE_var]

        # plot
        axs[0,0].pcolormesh(p_mod_m, vmin=vmin, vmax=vmax)
        axs[1,0].pcolormesh(p_mod_d, vmin=vmin, vmax=vmax)
        axs[0,1].pcolormesh(c_mod_m, vmin=vmin, vmax=vmax)
        axs[1,1].pcolormesh(c_mod_d, vmin=vmin, vmax=vmax)
        axs[0,2].pcolormesh(sat_m, vmin=vmin, vmax=vmax)
        axs[1,2].pcolormesh(sat_d, vmin=vmin, vmax=vmax)

        plt.savefig(f"FIGS/CO9_CO7_CMEMS_Satellite_{KE_var}.png", dpi=600)

    def plot_KE_master(self, vmin=0, vmax=0.06):
        """
        plot KE for satellite data and two models
        a 3x3 panel plot with:
           - row 0 : monthly MKE for all cases
           - row 1 : monthly EKE for all cases
           - row 2 : daily EKE for models
        """

        # get domain cfg
        domcfg = xr.open_dataset(cfg.dn_dom + cfg.grid_nc,
                                 chunks="auto").squeeze()

        # set projection
        mid_lat = np.mean([44,domcfg.nav_lat.max().values])
        mid_lon = -2.6
        proj=ccrs.EquidistantConic(central_latitude=mid_lat,
          standard_parallels=(44,domcfg.nav_lat.max().values),
          central_longitude=mid_lon) 
        plt_proj=ccrs.PlateCarree()
        proj_dict = {"projection": proj}

        # initialise figure
        fig, axs = plt.subplots(3, 3, figsize=(5.5,5.0), subplot_kw=proj_dict)
        plt.subplots_adjust(left=0.08, right=0.87, top=0.95, bottom=0.10,
                            hspace=0.01, wspace=0.04)

        # get primary model data
        path = f"{cfg.dn_out}/satellite/"
        p_mod_d = xr.open_dataset(path + "2004_2014_25hourm_KE.nc",
                                  chunks="auto")
        p_mod_m = xr.open_dataset(path + "2004_2014_monthly_KE.nc",
                                  chunks="auto")

        # get satellite data
        sat0_d = xr.open_dataset(
                      path + "2004-01-01_2014-01-01_1D_CMEMS_satellite_KE.nc",
                                  chunks="auto")
        sat0_m = xr.open_dataset(
                      path + "2004-01-01_2014-01-01_1MS_CMEMS_satellite_KE.nc",
                                  chunks="auto")

        # get comparison model data
        path = cfg.comp_case["proc_data"] + "/satellite/"
        c_mod_d = xr.open_dataset(path + "2004_2014_25hourm_KE.nc",
                                  chunks="auto")
        c_mod_m = xr.open_dataset(path + "2004_2014_monthly_KE.nc",
                                  chunks="auto")


        # plot
        def render_ke(ax, ds, KE_type):
        
            p = ax.pcolormesh(domcfg.nav_lon, domcfg.nav_lat, ds[KE_type],
                          vmin=vmin, vmax=vmax, transform=plt_proj,
                          cmap=plt.cm.Oranges)

            return p

        render_ke(axs[0,0], p_mod_m, "MKE")
        render_ke(axs[0,1], p_mod_m, "EKE")
        render_ke(axs[0,2], p_mod_d, "EKE")

        render_ke(axs[1,0], c_mod_m, "MKE")
        render_ke(axs[1,1], c_mod_m, "EKE")
        render_ke(axs[1,2], c_mod_d, "EKE")

        render_ke(axs[2,0], sat0_m, "MKE")
        p = render_ke(axs[2,1], sat0_m, "EKE")
        render_ke(axs[2,2], sat0_d, "EKE")

        for ax in axs[:,1:].flatten():
            ax.set_yticklabels([])
        for ax in axs[:-1].flatten():
            ax.set_xticklabels([])
        for ax in axs[:,0]:
            ax.set_ylabel("Latitude")
        for ax in axs[-1]:
            ax.set_xlabel("Longitude")

        #for ax in axs.flatten():
        #    ax.set_aspect("equal")

        axs[0,0].text(0.5, 1.05, "MKE",
                     rotation=0, transform=axs[0,0].transAxes,
                     va="bottom", ha="center")
        axs[0,1].text(0.5, 1.05, "EKE monthly data",
                     rotation=0, transform=axs[0,1].transAxes,
                     va="bottom", ha="center")
        axs[0,2].text(0.5, 1.05, "EKE daily data",
                     rotation=0, transform=axs[0,2].transAxes,
                     va="bottom", ha="center")

        for i in range(3):
            axs[0,i].text(0.05,0.95, "CO9", ha="left", va="top",
                          transform=axs[0,i].transAxes)
            axs[1,i].text(0.05,0.95, "CO7", ha="left", va="top",
                          transform=axs[1,i].transAxes)
            axs[2,i].text(0.05,0.95, "Satellite L4", ha="left", va="top",
                          transform=axs[2,i].transAxes)

        # colorbar
        pos0 = axs[0,-1].get_position()
        pos1 = axs[1,-1].get_position()
        pos2 = axs[2,-1].get_position()

        cbar_ax = fig.add_axes([0.88, 0.2, 
                                0.02, 0.6])

        cbar = fig.colorbar(p, cax=cbar_ax, orientation='vertical')
        cbar.ax.text(4.5, 0.5, r"Kinetic Energy (m$^2$ s$^{-1}$)",
                     rotation=90, transform=cbar.ax.transAxes,
                     va='center', ha='left')

        for i, ax in enumerate(axs.flatten()):

            if i in [0,1,2,3,4,5]:
                b_label=False
            else:
                b_label=True

            if i in [1,2,4,5,7,8]:
                l_label=False
            else:
                l_label=True

            self.add_land_and_gridlines(ax, domcfg, plt_proj, t_label=False,
                                        r_label=False, b_label=b_label,
                                        l_label=l_label)

        plt.savefig(f"FIGS/CO9_CO7_CMEMS_Satellite_KE.png", dpi=600)

    def add_land_and_gridlines(self, ax, domcfg, plt_proj,
                               l_label=True, r_label=True,
                               t_label=True, b_label=True):
        # add land mask
        landmask = xr.where(domcfg.bottom_level == 0, 1, np.nan)
        print (landmask)
        ax.contourf(domcfg.nav_lon, domcfg.nav_lat, landmask, 
                   colors=[plt.cm.Greys(0.2)],
                   transform=plt_proj)

        # format gridlines
        lon_grid = [-20,-10, 0, 10]
        lat_grid = [45, 50, 55, 60, 65]
        gl = ax.gridlines(draw_labels=True, xlocs=lon_grid, ylocs=lat_grid,
                         color='k', alpha=0.20)
        gl.xpadding = 2
        gl.ypadding = 2
        gl.xlabel_style = {'size': 8}
        gl.ylabel_style = {'size': 8}
        
        gl.bottom_labels = b_label
        gl.left_labels = l_label
        gl.right_labels = r_label
        gl.top_labels = t_label

        plt.draw()
    
    

    def plot_KE_product_comparison(self, vmin=0, vmax=0.06):
        """ side comparison of CMEMS and AVISO_SWOT data """
        # initialise figure
        fig, axs = plt.subplots(4, 3, figsize=(5.5,6.5))
        plt.subplots_adjust(left=0.1, right=0.88, top=0.98, bottom=0.1,
                            hspace=0.02, wspace=0.02)

        # get primary model data
        path = f"{cfg.dn_out}/satellite/"
        p_mod_d = xr.open_dataset(path + "2004_2014_25hourm_KE.nc",
                                  chunks="auto")
        p_mod_m = xr.open_dataset(path + "2004_2014_monthly_KE.nc",
                                  chunks="auto")

        # get satellite data
        #sat0_d = xr.open_dataset(
        #              path + "2004-01-01_2014-01-01_1D_CMEMS_satellite_KE.nc",
        #                          chunks="auto")
        sat0_m = xr.open_dataset(
                      path + "2004-01-01_2014-01-01_1MS_CMEMS_satellite_KE.nc",
                                  chunks="auto")

        # get swot data
        sat1_m = xr.open_dataset(
                path + "2023-03-01_2025-01-01_1MS_AVISO_SWOT_satellite_KE.nc",
                                  chunks="auto")
        sat1_d = xr.open_dataset(
                path + "2023-03-28_2025-01-11_1D_AVISO_SWOT_satellite_KE.nc",
                                  chunks="auto")

        # get comparison model data
        path = cfg.comp_case["proc_data"] + "/satellite/"
        c_mod_d = xr.open_dataset(path + "2004_2014_25hourm_KE.nc",
                                  chunks="auto")
        c_mod_m = xr.open_dataset(path + "2004_2014_monthly_KE.nc",
                                  chunks="auto")

        # get domain cfg
        domcfg = xr.open_dataset(cfg.dn_dom + cfg.grid_nc, chunks="auto")

        # plot
        def render_ke(ax, ds, KE_type):
        
            ax.pcolormesh(domcfg.nav_lon, domcfg.nav_lat, ds[KE_type],
                          vmin=vmin, vmax=vmax)

        render_ke(axs[0,0], p_mod_m, "MKE")
        render_ke(axs[0,1], p_mod_m, "EKE")
        render_ke(axs[0,2], p_mod_d, "EKE")

        render_ke(axs[1,0], c_mod_m, "MKE")
        render_ke(axs[1,1], c_mod_m, "EKE")
        render_ke(axs[1,2], c_mod_d, "EKE")

        render_ke(axs[2,0], sat0_m, "MKE")
        render_ke(axs[2,1], sat0_m, "EKE")
        #render_ke(axs[2,2], sat0_d, "EKE")

        render_ke(axs[3,0], sat1_m, "MKE")
        render_ke(axs[3,1], sat1_m, "EKE")
        render_ke(axs[3,2], sat1_d, "EKE")

        for ax in axs[:,1:].flatten():
            ax.set_yticklabels([])
        for ax in axs[:-1].flatten():
            ax.set_xticklabels([])
        for ax in axs[:,0]:
            ax.set_ylabel("Latitude")
        for ax in axs[-1]:
            ax.set_xlabel("Longitude")

        # colorbar
        #pos = axs.get_position()
        #cbar_ax = fig.add_axes([pos.x0, 0.12, 
        #                        pos.x1 - pos.x0, 0.02])
        #cbar = fig.colorbar(p, cax=cbar_ax, orientation='horizontal')
        #cbar.ax.text(0.5, -2.8, r"Temperature ($^{\circ}$C)", fontsize=8,
        #             rotation=0, transform=cbar.ax.transAxes,
        #             va='top', ha='center')
        
        #plt.show()
        plt.savefig(f"FIGS/AVISO_SWOT_CMEMS_Satellite_KE.png", dpi=600)

    def PDF_KE(self):
        """ probability density function of KE """

        fig, axs = plt.subplots(2,1, figsize=(6.5,3))

        # get primary model data
        path = f"{cfg.dn_out}/satellite/"
        p_mod_d = xr.open_dataset(path + "2004_2014_25hourm_KE.nc",
                                  chunks="auto")
        p_mod_m = xr.open_dataset(path + "2004_2014_monthly_KE.nc",
                                  chunks="auto")

        # get satellite data
        #sat_d = xr.open_dataset(path + "2004_2014_1D_satellite_KE.nc",
        #                          chunks="auto")
        sat_m = xr.open_dataset(path + "2004_2014_1MS_satellite_KE.nc",
                                  chunks="auto")

        # get comparison model data
        path = cfg.comp_case["proc_data"] + "/satellite/"
        c_mod_d = xr.open_dataset(path + "2004_2014_25hourm_KE.nc",
                                  chunks="auto")
        c_mod_m = xr.open_dataset(path + "2004_2014_monthly_KE.nc",
                                  chunks="auto")

        bins = np.logspace(-4,-2,10)
        p_mod_m["MKE"] = xr.where(p_mod_m.MKE > 1e-6, p_mod_m.MKE, np.nan)
        c_mod_m["MKE"] = xr.where(c_mod_m.MKE > 1e-6, c_mod_m.MKE, np.nan)
        sat_m["MKE"] = xr.where(sat_m.MKE > 1e-6, sat_m.MKE, np.nan)
        p_mod_m["EKE"] = xr.where(p_mod_m.EKE > 1e-6, p_mod_m.EKE, np.nan)
        c_mod_m["EKE"] = xr.where(c_mod_m.EKE > 1e-6, c_mod_m.EKE, np.nan)
        sat_m["EKE"] = xr.where(sat_m.EKE > 1e-6, sat_m.EKE, np.nan)

        axs[0].hist(p_mod_m.MKE.stack(z=["x","y"]), bins, density=True, alpha=0.4)
        axs[0].hist(c_mod_m.MKE.stack(z=["x","y"]), bins, density=True, alpha=0.4)
        axs[0].hist(sat_m.MKE.stack(z=["x","y"]), bins, density=True, alpha=0.4)

        axs[1].hist(p_mod_m.EKE.stack(z=["x","y"]), bins, density=True, alpha=0.4)
        axs[1].hist(c_mod_m.EKE.stack(z=["x","y"]), bins, density=True, alpha=0.4)
        axs[1].hist(sat_m.EKE.stack(z=["x","y"]), bins, density=True, alpha=0.4)

        for ax in axs:
            ax.set_xscale("log")

        plt.show()

def get_eof(ds, dn_out, fn, t_mode=False):
    """ calculate eof of surface data """

    # initiate eof model
    # Note: use_coslat should be used to weight but latitude, but xeof cannot
    # handle 2d latitude variable - it searches for coordinate dimensions
    #print (ds)
    #dsnan = np.isnan(ds)
    #for i in range(120):
    #    plt.pcolor(ds.isel(time=i))
    #    plt.colorbar()
    #    plt.show()
    #print (dsnan)
    #print (sdhfkj)
    
    model = xe.single.EOF(n_modes=5, standardize=False)

    # calculate eof
    ds = ds.assign_coords({"x":np.arange(ds.sizes["x"]),
                           "y":np.arange(ds.sizes["y"])})

    if t_mode:
        model.fit(ds, dim=("x","y"))
    else:
        model.fit(ds, dim=("time"))

    # save components to netcdf
    components = model.components(normalized=False)
    del components.attrs["solver_kwargs"]  # attr causes error
    components.to_netcdf(f"{dn_out}/satellite/{fn}_eof_components.nc")

    # save scores to netcdf
    scores = model.scores(normalized=False)
    del scores.attrs["solver_kwargs"]  # attr causes error
    scores.to_netcdf(f"{dn_out}/satellite/{fn}_eof_scores.nc")

    # save explained variance to netcdf
    var_exp = model.explained_variance_ratio()
    del var_exp.attrs["solver_kwargs"]  # attr causes error
    var_exp.to_netcdf(f"{dn_out}/satellite/{fn}_eof_var_explained_ratio.nc")

if __name__ == "__main__":

    def get_co9_gridded_satellite_data(var="sst", time_step="1MS"):
        mod_cfg_fn = cfg.dn_dom + cfg.grid_nc
        sat = satellite(mod_cfg_fn)
        sat.get_cmems(var=var)
        if time_step == "1MS":
            sat.monthly_mean()
        with ProgressBar():
            sat.ds = sat.ds.load()
        sat.ds = sat.interpolate_to_model(sat.ds)
        sat.save_ds(f"CMEMS_L4_satellite_{var}_{time_step}_gridded_to_{cfg.case}.nc")
    def get_co9_raw_satellite_data(var="sst", time_step="1MS"):
        mod_cfg_fn = cfg.dn_dom + cfg.grid_nc
        sat = satellite(mod_cfg_fn)
        sat.get_cmems(var=var)
        if time_step == "1MS":
            sat.monthly_mean()
        with ProgressBar():
            sat.ds = sat.ds.load()
        sat.save_ds(f"CMEMS_L4_satellite_{var}_{time_step}_raw.nc")

    def grid_co9_raw_satellite_data(var, time_step):
        mod_cfg_fn = cfg.dn_dom + cfg.grid_nc
        sat = satellite(mod_cfg_fn)
        fn_name=f"CMEMS_L4_satellite_{var}_{time_step}_raw.nc"
        fn = f"{cfg.dn_out}/satellite/{fn_name}"
        sat.ds = xr.open_dataarray(fn, chunks="auto")
        print (sat.ds)
        with ProgressBar():
            sat.ds = sat.ds.load()

        sat.var_str = var
        sat.ds = sat.interpolate_to_model(sat.ds)
        sat.save_ds(f"CMEMS_L4_satellite_{var}_{time_step}_gridded_to_{cfg.case}.nc")

    def get_co9_gridded_satellite_data_inc_swot(time_step="1MS"):
        mod_cfg_fn = cfg.dn_dom + cfg.grid_nc
        sat = satellite(mod_cfg_fn)
        sat.get_l4_SWOT()
        if time_step == "1MS":
            sat.monthly_mean()
        with ProgressBar():
            sat.ds = sat.ds.load()
        ds_set= []
        for var in ["ugos","vgos"]:
            sat.var_str = var
            ds_set.append(sat.interpolate_to_model(sat.ds[var]))
        sat.ds = xr.merge(ds_set)
        sat.save_ds(f"dt_allsat_phy_l4_{time_step}_gridded_to_{cfg.case}.nc")

    def calculate_satellite_eof(var_nam="ssh", fn_nam="ssh"):
        
        path = f"{cfg.dn_out}/satellite/"
        fn = f"{path}/CMEMS_L4_satellite_{var}_gridded_to_{cfg.case}.nc"
        sat_proc = xr.open_dataarray(fn, chunks=-1)
        sat_proc["time"] = sat_proc.time.astype("datetime64[M]")

        # remove deep water
        cfg_fn = '/gws/nopw/j04/jmmp/public/AMM15/DOMAIN_CFG/GEG_SF12.nc'
        domcfg = xr.open_dataset(cfg_fn).squeeze()

        sat_proc = sat_proc.where(domcfg.bathy < 200)

        # remove shallow atlantic areas
        sat_proc = sat_proc.where(sat_proc.longitude>-12.5)
        sat_proc = sat_proc.where((sat_proc.longitude>-4.5) | 
                                  (sat_proc.latitude<60))

        get_eof(sat_proc, cfg.dn_out, f"CMEMS_L4_satellite_{fn_nam}") 


    def calculate_primary_model_eof(var_nam="sossheig", fn_nam="ssh"):

        # get model and remove surface loading 
        fn = cfg.dn_dat
        mod = model_surface(fn, src_t_coord="time_counter",
                                src_x_coord="x_grid_T",
                                src_y_coord="y_grid_T")
        mod.cfg_fn = '/gws/nopw/j04/jmmp/public/AMM15/DOMAIN_CFG/GEG_SF12.nc'
        mod.get_time_mean_var(var_nam=var_nam)
        mod.map_lat_lon_names("nav_lon_grid_T", "nav_lat_grid_T")
        mod.remove_inverse_barometer()

        # retrieve dataset
        mod_proc = mod.ds

        # remove deep water
        domcfg = mod.get_domain_cfg()
        mod_proc = mod_proc.where(domcfg.bathy < 200)
        
        # remove shallow atlantic areas
        mod_proc = mod_proc.where(mod_proc.longitude>-12.5)
        mod_proc = mod_proc.where((mod_proc.longitude>-4.5) | 
                                  (mod_proc.latitude<60))

        # get eof of ssh
        get_eof(mod_proc, cfg.dn_out, f"CO9_{fn_nam}")


    def calculate_comparison_model_eof(var_nam="sossheig", fn_nam="ssh"):

        # get model and remove surface loading 
        fn_raw = cfg.comp_case["raw_data"]
        fn_proc = cfg.comp_case["proc_data"]
        mod = model_surface(fn_raw, fn_proc, src_t_coord="time_counter")
        mod.cfg_fn = cfg.dn_dom + cfg.comp_case["grid"]
        mod.get_time_mean_var(var_nam=var_nam)
        mod.map_lat_lon_names("nav_lon", "nav_lat")
        mod.remove_inverse_barometer("era_interim")

        # retrieve dataset
        mod_proc = mod.ds

        # remove deep water
        domcfg = mod.get_domain_cfg(rename="hbatt")
        mod_proc = mod_proc.where(domcfg.bathy < 200)

        # remove shallow atlantic areas
        mod_proc = mod_proc.where(mod_proc.longitude>-12.5)
        mod_proc = mod_proc.where((mod_proc.longitude>-4.5) | 
                                  (mod_proc.latitude<60))

        #plt.pcolormesh(m.nav_lon, m.nav_lat,m, cmap=plt.cm.binary)
        #plt.pcolormesh(m.nav_lon, m.nav_lat,m_cut)
        #plt.pcolormesh(m.nav_lon, m.nav_lat,m_cut_f)
        #plt.axvline(-12.5)
        #plt.show()
        #mod_proc.squeeze().plot()
        
        # get eof of ssh
        get_eof(mod_proc, cfg.comp_case["proc_data"],
                cfg.comp_case["case"] + "_{var_nam}")

    def plot_eof():
        splot = satellite_plot()
        splot.plot_eof_validation("CO9","co7", "CMEMS_L4_satellite", "sst")

    def save_var_primary_model(var_nam, grid):
        fn_raw = cfg.dn_dat
        fn_proc = cfg.dn_out
        mod = model_surface(fn_raw, fn_proc, src_t_coord="time_counter",
                                             src_z_coord="depth" + grid.lower())
        mod.get_time_mean_var(var_nam=var_nam, freq=None, grid=grid, save=True)

    def save_var_comparison_model(var_nam, grid):
        fn_raw = cfg.comp_case["raw_data"]
        fn_proc = cfg.comp_case["proc_data"]
        mod = model_surface(fn_raw, fn_proc, src_t_coord="time_counter",
                            src_z_coord="depth" + grid.lower())
        mod.get_time_mean_var(var_nam=var_nam, freq=None, grid=grid, save=True)

    def get_KE(fn_raw, fn_proc, freq="monthly"):
        mod = model_surface(fn_raw, fn_proc, src_t_coord="time")

        u = xr.open_dataarray(fn_proc + f"{cfg.y0}_{cfg.y1}_{freq}_vozocrtx.nc",
                chunks="auto")
        v = xr.open_dataarray(fn_proc + f"{cfg.y0}_{cfg.y1}_{freq}_vomecrty.nc",
                chunks="auto")
        KE = mod.get_KE(u, v)
        KE.to_netcdf(fn_proc + f"satellite/{cfg.y0}_{cfg.y1}_{freq}_KE.nc")

    def get_satellite_KE():
        mod_cfg_fn = cfg.dn_dom + cfg.grid_nc
        sat = satellite(mod_cfg_fn)
        sat.get_KE(save=True, time_step="1D", src="CMEMS")
        #sat.get_KE(save=True, time_step="1MS", src="CMEMS")
        

    #get_satellite_KE()
    #sp = satellite_plot()
    #sp.PDF_KE()
    #sp.plot_KE_master()
    #get_KE(cfg.comp_case["raw_data"], cfg.comp_case["proc_data"], freq="25hourm")
    #get_KE(cfg.dn_dat, cfg.dn_out, freq="25hourm")
    #get_KE(cfg.comp_case["raw_data"], cfg.comp_case["proc_data"])
    #get_KE(cfg.dn_dat, cfg.dn_out)
    #save_var_primary_model("vomecrty", "V")
    #save_var_primary_model("vozocrtx", "U")
    #save_var_comparison_model("vomecrty", "V")
    #save_var_comparison_model("vozocrtx", "U")
    
    plot_eof()
    #calculate_primary_model_eof("vomecrx")
    #calculate_comparison_model_eof("votemper")
    #calculate_comparison_model_eof("vosaline")
    #calculate_satellite_eof()
    #get_co9_gridded_satellite_data_inc_swot(time_step="1MS")
    #get_co9_gridded_satellite_data_inc_swot(time_step="1D")
    #get_co9_gridded_satellite_data(var="uvel", time_step="1MS")
    #get_co9_gridded_satellite_data(var="vvel", time_step="1MS")
    #get_co9_gridded_satellite_data(var="uvel", time_step="1D")
    #get_co9_gridded_satellite_data(var="vvel", time_step="1D")
    #grid_co9_raw_satellite_data(var="vvel", time_step="1D")
    #get_co9_raw_satellite_data(var="vvel", time_step="1D")
