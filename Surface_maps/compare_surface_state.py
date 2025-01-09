from PythonEnvCfg.config import config, bounds
cfg = config() # initialise variables in python

import coast
import xarray as xr
import os
import copernicusmarine
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import numpy as np

class extract_surface(object):
    def __init__(self):

        # paths
        self.fn_dom = config.dn_dom + config.grid_nc
        self.fn_dat = config.dn_out + "profiles/gridded*.nc"
        self.fn_out = config.dn_out + 'surface_maps/'

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
           data_request["fn"] = "cmems_obs-sst_atl_phy_nrt_l3s_P1D-m",
           data_request["variables"] = ["sea_surface_temperature"]

        if var == "ssh":
           data_request["fn"] = "cmems_obs-sl_eur_phy-ssh_my_allsat-l4-duacs-0.0625deg_P1D"
           data_request["variables"] = ["adt"]

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
        )
        #self.ds = self.ds.transpose("time","longitude","latitude")

    def monthly_mean(self):
        """ average over month """

        self.ds = self.ds.resample(time='1ME').mean().load()

    def interpolate_to_model(self, cfg_fn):
        """ interpolate lat-lon to horizontal grid """

        self.ds = self.ds.isel(time=slice(None,5))
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
            values = (ds_t.to_dataarray().values.flatten())
            
            n_grid.append(
             griddata(points, values, target, method="nearest")[:,:,np.newaxis])

        n_grid_all = np.concatenate(n_grid, axis=2)

        self.ds = xr.DataArray(
                             data=n_grid_all,
                             dims=["x","y","time"],
                             coords={"longitude": (["x","y"],tgt_lon.values),
                                     "latitude": (["x","y"],tgt_lat.values),
                                     "time": self.ds.time},
                             name="adt").to_dataset()


    def quick_compare(self):
        """ quick plot """

        snapshot_obs = self.ds.isel(time=0).adt
        snapshot_co9 = self.ds.isel(time=0)

        plt.pcolor(snapshot)
        plt.show()

class model_surface(object):

    def __init__(self, fn_path):

        self.ds = xr.open_mfdataset(fn_path + "*.nc.ppc3", chunks=-1).zos

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

if __name__ == "__main__":

    # get satellite
    sat = satellite()
    sat.get_cmems()
    sat.monthly_mean()
    cfg_fn = '/gws/nopw/j04/jmmp/public/AMM15/DOMAIN_CFG/GEG_SF12.nc'
    sat.interpolate_to_model(cfg_fn)

    # get model
    fn = "/gws/nopw/j04/jmmp/jmmp_collab/AMM15/OUTPUTS/P1.5c/MONTHLY/"
    mod = model_surface(fn)

    # model
    splot = satellite_plot()
    splot.plot_model_and_satellite_snapshot_ssh(mod.ds, sat.ds.adt)
