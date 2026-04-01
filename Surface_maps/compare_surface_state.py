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

    def monthly_mean(self):
        """ average over month """

        self.ds = self.ds.resample(time='1MS').mean().load()

    def interpolate_to_model(self, cfg_fn):
        """ interpolate lat-lon to horizontal grid """

        domcfg = xr.open_dataset(cfg_fn)

        tgt_lon =  domcfg.nav_lon
        tgt_lat =  domcfg.nav_lat
        
        target = (tgt_lon, tgt_lat)

        src_lon = self.ds.longitude.values
        src_lat = self.ds.latitude.values

        src_mlon, src_mlat = np.meshgrid(src_lon, src_lat)

        points = (src_mlon.flatten(), src_mlat.flatten())

        n_grid = []
        for time, ds_t in self.ds.groupby("time"):
            print (time)
            values = ds_t.values.flatten()
            values_masked = (np.nan_to_num(values))
            
            n_grid.append(
             griddata(points, values_masked, target,
                      method="cubic")[:,:,np.newaxis])

        n_grid_all = np.concatenate(n_grid, axis=2)

        self.ds = xr.DataArray(
                             data=n_grid_all,
                             dims=["y","x","time"],
                             coords={"longitude": (["y","x"],tgt_lon.values),
                                     "latitude": (["y","x"],tgt_lat.values),
                                     "time": self.ds.time},
                             name=self.var_str)#.to_dataset()

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

    def get_KE(self, save=False):
        """ get mean and eddy kinetic energy of surface currents """

        path = f"{cfg.dn_out}/satellite/"
        fn = path + "CMEMS_L4_satellite_uvel_gridded_to_P2.0.nc"
        u = xr.open_dataarray(fn)
        fn = path + "CMEMS_L4_satellite_vvel_gridded_to_P2.0.nc"
        v = xr.open_dataarray(fn)

        u_bar = u.mean("time")
        v_bar = v.mean("time")

        u_prime = u_bar - u
        v_prime = v_bar - v

        MKE = 0.5 * (u_bar**2 + v_bar**2)
        EKE = 0.5 * ((u_prime**2).mean("time") +
                     (v_prime**2).mean("time"))

        MKE.name = "MKE"
        EKE.name = "EKE"

        KE = xr.merge([MKE,EKE])

        if save: 
            # needs reduced hardcoding
            fn = f"{cfg.dn_out}/satellite/2004_2014_satellite_KE.nc"
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

    def plot_eof_validation(self, mod0, mod1, sat, var):
        """ plot eof breakdown of model versus obs """

        # initialise figure
        fig = plt.figure(figsize=(5.5,6.5))

        # initialise gridspec
        gs0 = gridspec.GridSpec(ncols=2, nrows=1)
        gs1 = gridspec.GridSpec(ncols=2, nrows=3)

        ## set frame bounds
        gs0.update(top=0.96, bottom=0.78, left=0.1, wspace=0.1, hspace=0.12,
                   right=0.85)
        gs1.update(top=0.7, bottom=0.08, left=0.1, wspace=0.1, hspace=0.08,
                   right=0.85)

        # set projection
        proj=ccrs.AlbersEqualArea()
        proj=ccrs.PlateCarree()

        # assign axes to lists
        axs0 = []
        for i in range(2):
            axs0.append(fig.add_subplot(gs0[i]))
        axs1 = []
        for i in range(6):
            axs1.append(fig.add_subplot(gs1[i], projection=proj))

        plt_proj=ccrs.PlateCarree()
        proj_dict = {"projection": proj}

        # get data
        path = f"{cfg.dn_out}/satellite/"
        eof_type = ""
        mod0_scores = xr.open_dataarray(f"{path}{mod0}_{var}_eof{eof_type}_scores.nc")
        mod0_comp = xr.open_dataarray(f"{path}{mod0}_{var}_eof{eof_type}_components.nc")

        sat_scores = xr.open_dataarray(f"{path}{sat}_{var}_eof{eof_type}_scores.nc")
        sat_comp = xr.open_dataarray(f"{path}{sat}_{var}_eof{eof_type}_components.nc")

        #path = cfg.comp_case["proc_data"] + "/satellite/"
        mod1_scores = xr.open_dataarray(f"{path}{mod0}_{var}_eof{eof_type}_scores.nc")
        mod1_comp = xr.open_dataarray(f"{path}{mod0}_{var}_eof{eof_type}_components.nc")

        # get p-value
        mod0_pval = xr.corr(mod0_scores, sat_scores, dim="time")
        mod1_pval = xr.corr(mod1_scores, sat_scores, dim="time")

        def render(axs, comp, scores, i, label):
            axs0[0].plot(scores.time, scores.sel(mode=1), label=label,
                          lw=0.8)
            axs0[1].plot(scores.time, scores.sel(mode=2), label=label,
                          lw=0.8)

            axs1[i*2].pcolormesh(comp.longitude, comp.latitude,
                                comp.sel(mode=1).squeeze())
            p = axs1[i*2+1].pcolormesh(comp.longitude, comp.latitude,
                                comp.sel(mode=2).squeeze())
            #for ax in axs1[i*2:i*2+2]:
            #    ax.pcolormesh(comp.longitude, comp.latitude,
            #        xr.where(comp.isel(mode=1).squeeze() == np.nan, 1, np.nan))
            plt.colorbar(p)


            # set extent
            #lon0 = comp.longitude.isel(x=0, y=0).values
            lon0 = -15
            lon1 = 9.8
            #axs[i,0].set_extent([lon0, lon1, 46, 62])
            #axs[i,1].set_extent([lon0, lon1, 46, 62])
            axs1[i].set_xlim(lon0, lon1)
            axs1[i].set_ylim(46, 62)
            axs1[i+3].set_xlim(lon0, lon1)
            axs1[i+3].set_ylim(46, 62)

        render(axs1, mod0_comp, mod0_scores, 0, "CO9")
        render(axs1, mod1_comp, mod1_scores, 1, "CO7")
        print (sat_comp)
        #sat_comp = sat_comp.rename({"nav_lon":"longitude",
        #                            "nav_lat":"latitude"})
        render(axs1, sat_comp, sat_scores, 2, "Obs")

        axs0[1].legend(loc="upper left", bbox_to_anchor=(1.02,1),
                        bbox_transform=axs0[1].transAxes)

        # set timeseries lims
        for ax in axs0:
            ax.set_xlim(mod0_scores.time.min(), mod0_scores.time.max())

        # add p-vals
        p =  str(np.round(mod0_pval.sel(mode=1).data, 2))
        axs0[0].text(0.5, 0.95, "p = " + p,
                      ha="left", va="top", transform=axs0[0].transAxes)
        p =  str(np.round(mod1_pval.sel(mode=1).data, 2))
        axs0[0].text(0.75, 0.95, "p = " + p,
                      ha="left", va="top", transform=axs0[0].transAxes)
        p =  str(np.round(mod0_pval.sel(mode=2).data, 2))
        axs0[1].text(0.5, 0.95, "p = " + p,
                      ha="left", va="top", transform=axs0[1].transAxes)
        p =  str(np.round(mod1_pval.sel(mode=2).data, 2))
        axs0[1].text(0.75, 0.95, "p = " + p,
                      ha="left", va="top", transform=axs0[1].transAxes)

        # set labels
        axs0[0].set_title("mode 1")
        axs0[1].set_title("mode 2")

        axs1[0].text(0.05,0.95, "CO9", ha="left", va="top",
                      transform=axs1[0].transAxes)
        axs1[2].text(0.05,0.95, "CO7", ha="left", va="top",
                      transform=axs1[2].transAxes)
        axs1[4].text(0.05,0.95, "Obs", ha="left", va="top",
                      transform=axs1[4].transAxes)

        axs1[1].text(0.05,0.95, "CO9", ha="left", va="top",
                      transform=axs1[1].transAxes)
        axs1[3].text(0.05,0.95, "CO7", ha="left", va="top",
                      transform=axs1[3].transAxes)
        axs1[5].text(0.05,0.95, "Obs", ha="left", va="top",
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

        #cfg_fn = '/gws/nopw/j04/jmmp/public/AMM15/DOMAIN_CFG/GEG_SF12.nc'
       # domcfg = xr.open_dataset(cfg_fn).squeeze()
       # print (domcfg)
       # msk = xr.where(domcfg.top_level == 0, 0, np.nan)
       # print (msk)
        for ax in axs1:
            #ax.pcolormesh(domcfg.nav_lon, domcfg.nav_lat, msk, cmap="grey")
            #ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')
            land_50m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                        edgecolor='k',
                                        facecolor='grey')
            ax.add_feature(land_50m, zorder=100, lw=0.5)

            ax.set_xticks([-15, -10, -5, 0, 5, 10], crs=ccrs.PlateCarree())
            ax.set_yticks([50, 55, 60], crs=ccrs.PlateCarree())
            lon_formatter = LongitudeFormatter(zero_direction_label=True)
            lat_formatter = LatitudeFormatter()

        eof_str = eof_type.lstrip("_")
        print (eof_str)
        plt.savefig(f"FIGS/CO9_CO7_CMEMS_Satellite_{var}_{eof_str}_pca.png",
                     dpi=600)

    def plot_KE():
        """
        plot KE for satellite data and two models
        a 2x3 panel plot with:
           - row 1 : monthly EKE for all cases
           - row 2 : daily EKE for models
        """
        # initialise figure
        fig = plt.figure(figsize=(5.5,6.5))

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

    def get_co9_gridded_satellite_data(var="sst"):
        sat = satellite()
        sat.get_cmems(var=var)
        sat.monthly_mean()
        cfg_fn = '/gws/nopw/j04/jmmp/public/AMM15/DOMAIN_CFG/GEG_SF12.nc'
        sat.interpolate_to_model(cfg_fn)
        sat.save_ds(f"CMEMS_L4_satellite_{var}_gridded_to_{cfg.case}.nc")

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
        sat = satellite()
        sat.get_KE(save=True)
        

    get_satellite_KE()
    #get_KE(cfg.comp_case["raw_data"], cfg.comp_case["proc_data"], freq="25hourm")
    #get_KE(cfg.dn_dat, cfg.dn_out, freq="25hourm")
    #get_KE(cfg.comp_case["raw_data"], cfg.comp_case["proc_data"])
    #get_KE(cfg.dn_dat, cfg.dn_out)
    #save_var_primary_model("vomecrty", "V")
    #save_var_primary_model("vozocrtx", "U")
    #save_var_comparison_model("vomecrty", "V")
    #save_var_comparison_model("vozocrtx", "U")
    
    #plot_eof()
    #calculate_primary_model_eof("vomecrx")
    #calculate_comparison_model_eof("votemper")
    #calculate_comparison_model_eof("vosaline")
    #calculate_satellite_eof()
    #get_co9_gridded_satellite_data("uvel")
    #get_co9_gridded_satellite_data("vvel")
 
