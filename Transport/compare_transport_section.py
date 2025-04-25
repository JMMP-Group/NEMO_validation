from PythonEnvCfg.config import config
import xarray as xr
import pandas as pd
import numpy as np
import coast
from scipy.interpolate import NearestNDInterpolator
import matplotlib.pyplot as plt

cfg = config() # initialise variables in python

import time
class Ellet(object):
    def preprocess_annual_cruise_data(self):
        def open_ds(fn, path):
        
            df = pd.read_csv(path + fn, engine="c")
            ds = xr.Dataset.from_dataframe(df)
            
            return ds
        
        #path = "../../analysis_eel_data/data/raw/csv_datagridded/"
        path = "../../analysis_eel_data/data/raw/csv_ctdgrid/"
        data = open_ds("EELCTDandLADCP_3Dfield.csv", path)
    
        pos = open_ds("EELCTDandLADCP_refpos.csv", path)
        time = open_ds("EELCTDandLADCP_refdate.csv", path)

        
        # unstack dist, year and depth
        data["Refdist"] = data.Refdist.astype("int")
        data = data.set_coords(["Refdist","Year","Depth"])
        data = data.set_xindex(["Refdist","Year","Depth"])
        data = data.unstack("index")
        
        # set distance to index
        pos["Refdist"] = pos.Refdist.astype("int")
        pos = pos.set_coords("Refdist")
        pos = pos.swap_dims({"index":"Refdist"})
        pos = pos.drop_vars("index")
        
        # set time to index
        time = time.set_coords("Year")
        time = time.swap_dims({"index":"Year"})
        time = time.drop_vars("index")

        ds = xr.merge([data,pos,time])

        # set lat, lon as coords
        self.ds = ds.set_coords(["LonSta","LatSta"])
    
        # reduce to decade
        self.ds = self.ds.sel(Year=slice("2003","2014"))
        self.ds = self.ds.swap_dims({"Depth":"z_dim"})
        #self.ds["Depth"], _ = xr.broadcast(self.ds.Depth, self.ds.CruiseID)
        #print (self.ds)
        #print (ksdj)

        # reduce to rockall trough - two step to maintain dims
        #self.ds = self.ds.where((self.ds.LonSta < -9) & (self.ds.LonSta > -14),
        #                        drop=True)
        # using sel in place of where avoids uncesesary broadcasting of Month
        ind = self.ds.Refdist.where((self.ds.LonSta < -9) & 
                                    (self.ds.LonSta > -14),
                                    drop=True)
        self.ds = self.ds.sel(Refdist=ind)
        self.get_dz()
        self.get_dx()

        #self.ds = self.ds.stack(id_dim=["Year","Refdist"])

    def get_dz(self):
        """
        get integral distances
        """

        dz = self.ds.Depth.data[1:] - self.ds.Depth.data[:-1]
        dz_0 = [self.ds.Depth.data[0] + dz[0]]
        dz_end = [dz[-1]*2]
        dz = dz[1:]/2 + dz[:-1]/2
        self.ds['dz'] = xr.DataArray(np.concatenate((dz_0, dz, dz_end)),
                                     dims=("z_dim"))
    def get_dx(self):
        """
        get integral distances
        """

        dx = self.ds.Refdist.data[1:] - self.ds.Refdist.data[:-1]
        dx_0 = [dx[0]]
        dx_end = [dx[-1]]
        dx = dx[1:]/2 + dx[:-1]/2
        self.ds['dx'] = xr.DataArray(np.concatenate((dx_0, dx, dx_end)),
                                     dims=("Refdist"))
    def get_volume_transport(self):
        """
        get volume transport
        """

        ds = self.Ellet_profiles.dataset
        dims = ["id_dim","z_dim"]
        ds["volume_transport"] = (ds.ladcp_velocity * ds.dx * ds.dz).sum(dims)

    def get_dates(self):

        dates = []
        fn_dates = []
        for year, y_ds in self.ds.groupby("Year"):
            dates.append(str(year) + "-" + str(y_ds.Month.data[0]).zfill(2))
        print (dates)
        self.ds["time"] =  np.array(dates, dtype="datetime64")

        for year, y_ds in self.ds.groupby("Year"):
            fn_dates.append(str(year) + str(y_ds.Month.data[0]).zfill(2))
        
        return fn_dates

    def coast_formatting(self):
        """ format observations into COAsT profile object """

        # move data in profile object
        Ellet_cfg = "./Ellet_cfg.json"
        self.Ellet_profiles = coast.Profile(config=Ellet_cfg)
        self.Ellet_profiles.dataset = self.ds
        self.Ellet_profiles.apply_config_mappings()

    def save_processed_ellet(self):
        """ save processed ellet line data """

        #self.Ellet_profiles.dataset = self.Ellet_profiles.dataset.reset_index("id_dim")
        path = cfg.dn_out + "transport/obs_for_ellet_line.nc"
        self.Ellet_profiles.dataset.to_netcdf(path)

class ModelVels(object):

    def __init__(self, dates):
        self.dates = dates
        paths = [cfg.dn_dat + date + "01T0000Z_daily_grid_U.nc" for date in
                 dates]
        #paths_V = [cfg.dn_dat + date + "01T0000Z_daily_grid_V.nc" for date in
        #         dates]
        #print (paths_U)
        #nemo = coast.Gridded(paths, cfg.dn_dom + cfg.grid_nc,
        #                     config=cfg.fn_cfg_nemo, multiple=True)
        #nemo.dataset = self.restrict_lat_lon(nemo.dataset)
        #print (nemo.dataset.time)
        #print (skdjf)
#
#        #nemo.dataset = nemo.dataset.resample("time.year_month").mean()
#        print (nemo.dataset.time)
#        print (skdjf)
#        nemo.dataset = nemo.dataset.vozocrtx.groupby("year_month").mean()
#        print (nemo.dataset)
#        print (skdjf)
        #ds = xr.open_mfdataset(paths).vozocrtx.groupby("year_month").mean()
        #V_mean = xr.open_mfdataset(paths_V).vomecrty.groupby("year_month").mean()
        #print (U_mean.time_counter.values)

        #self.vels = {"U":xr.open_mfdataset(paths_U).vozocrtx,
        #             "V":xr.open_mfdataset(paths_V).vomecrty}

        
    
    def interpolate_vec_to_obs(self, obs, vec_str, var):

        # loop over years - accessing single month files
        interp_list = []
        for i, date in enumerate(self.dates):
            print (date)

            # access data
            path = cfg.dn_dat + date + f"01T0000Z_daily_grid_{vec_str}.nc"
            ds = xr.open_dataset(path, chunks="auto")[var]

            # reduce - time mean and lat-lon slice
            if vec_str=="U":
                ds = ds.isel(depthu=0)
                ds = ds.rename({"depthu":"depth"})
            else:
                ds = ds.isel(depthv=0)
                ds = ds.rename({"depthv":"depth"})
            ds = ds.mean("time_counter").load()
            ds = ds.where((ds.nav_lon < -8).compute() &
                          (ds.nav_lon > -14).compute() &
                          (ds.nav_lat > 56.5).compute() &
                          (ds.nav_lat < 58).compute(), drop=True)

            # get interpolation weights
            if i == 0:
                points = self.get_src_interp_points(ds)
                
            # get matching obs time
            obs_3d = obs.dataset.sel(time=int(date[:4]))

            # interpolate
            interpolated = self.interpolate_to_obs_loc(ds, obs_3d, points)
            ds_interp = xr.DataArray(data=interpolated[np.newaxis],
                                coords={"time":("time", [int(date)]),
                              "depth":(["z_dim","Refdist"],obs_3d.depth.data),
                              "Refdist":(["Refdist"],obs_3d.Refdist.data)},
                                       dims=("time","z_dim","Refdist"),
                                       name=vec_str)
            interp_list.append(ds_interp)

        obs_loc_vel = xr.concat(interp_list, dim="time")

        return obs_loc_vel

    def interpolate_vels_to_obs(self, obs):

        U = self.interpolate_vec_to_obs(obs, "U", "vozocrtx")
        V = self.interpolate_vec_to_obs(obs, "V", "vomecrty")

        vels = xr.merge([U,V])

        path = cfg.dn_out + "transport/model_vels_for_ellet_line.nc"
        vels.to_netcdf(path)

    def get_src_interp_points(self, src):

        src_lon_3d = np.broadcast_to(src.nav_lon.data, src.shape)
        src_lat_3d = np.broadcast_to(src.nav_lat.data, src.shape)
        src_dep_3d = np.transpose(
                     np.broadcast_to(src.depth.data, src.shape[::-1]))

        # flatten input
        src_lon = src_lon_3d.flatten()
        src_lat = src_lat_3d.flatten()
        src_dep = src_dep_3d.flatten()
        

        data_bool = ~np.isnan(src.values.flatten())
        src_lon = src_lon[data_bool]
        src_lat = src_lat[data_bool]
        src_dep = src_dep[data_bool]
        points = list(zip(src_dep, src_lat, src_lon))

        return points

    def interpolate_to_obs_loc(self, src, obs, points):
        """ interpolate to cruise stations """
        # need to make this 2d

        values = src.values.flatten()
        values = values[~np.isnan(values)]

        tgt_lon =  np.broadcast_to(obs.longitude,obs.ladcp_velocity.shape[::-1])
        tgt_lat =  np.broadcast_to(obs.latitude, obs.ladcp_velocity.shape[::-1])
        tgt_dep = obs.depth

        interp = NearestNDInterpolator(points, values)
        n_grid = interp(tgt_dep, tgt_lat, tgt_lon)
        return n_grid

    def restrict_lat_lon(self, ds, lat_bounds=[56.5,58], lon_bounds=[-8,-14]):
        ds = ds.where((ds.longitude < lon_bounds[0]).compute() &
                      (ds.longitude > lon_bounds[1]).compute() &
                      (ds.latitude > lat_bounds[0]).compute() &
                      (ds.latitude < lat_bounds[1]).compute(), drop=True)
        return ds

    def interpolate_vels_to_obs_coast(self, obs):

        print ("pre grid")
        model_path = cfg.dn_dat + f"*01T0000Z_daily_grid_U.nc"
        nemo = coast.Gridded(model_path, cfg.dn_dom + cfg.grid_nc,
                             config=cfg.fn_cfg_nemo, multiple=True)

        # limit lat lon for efficiency
        nemo.dataset = self.restrict_lat_lon(nemo.dataset)

        nemo.dataset = nemo.dataset.resample(time="ME").mean()

        nemo.dataset["landmask"] = nemo.dataset.bottom_level == 0
        nemo.dataset = nemo.dataset.rename({"depth_0": "depth"})
        nemo_profiles = obs.obs_operator(nemo)


def process_ellet_obs():
    ell = Ellet()
    ell.preprocess_annual_cruise_data()
    ell.coast_formatting()
    ell.get_volume_transport()
    ell.save_processed_ellet()
process_ellet_obs()

def get_model_on_ellet_locs():
    ell = Ellet()
    ell.preprocess_annual_cruise_data()
    ell.coast_formatting()
    
    m = ModelVels(ell.get_dates())
    m.interpolate_vels_to_obs(ell.Ellet_profiles)

#get_model_on_ellet_locs()

