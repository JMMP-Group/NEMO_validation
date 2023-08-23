### work in progress ###
### replication of  EN4_processing/GEN_MOD_Dave_example_profile_validation.py
### purpose > make gridded version of en4 data?

from config import config
from config import bounds
import sys

cfg = config() # initialise variables in python
bdy = bounds("AMM15")

from NEMO_validation._utils import landmask
import coast
import xarray as xr
import numpy as np
import pandas as pd
import os
from coast import crps_util as cu
import numpy as np

import time

class gridded_en4(object):

    def __init__(self):

        # user inputs
        args = sys.argv
        startyear=int(args[1])
        month=int(args[2])
        endyear=int(args[3])

        # File paths 
        self.fn_dom = cfg.dn_dom + cfg.grid_nc
        self.fn_dat = "%s%s%02d*T.nc"%(cfg.dn_dat, startyear, month) 
        self.fn_out = cfg.dn_out + 'surface_maps/'
        
        # Make out directory
        print(os.popen(f"mkdir -p {cfg.dn_out}").read())
        
        # Generated by pre_process_en4_monthly.py
        #fn_append = "_processed_%d%02d.nc"%(startyear,month)
        fn_append = "_processed_*.nc"
        self.fn_prof = cfg.dout_en4 + cfg.region + fn_append
        self.en4_path = cfg.dout_en4 + cfg.region 
        
        self.start_date = np.datetime64(str(startyear)+"-01-01")
        self.end_date = np.datetime64(str(endyear)+"-01-01")

        self.month_list = np.arange(self.start_date, self.end_date, 
                                    dtype="datetime64[M]")
        
    def model_data(self):
        # CREATE NEMO OBJECT and read in NEMO data.
        nemo = coast.Gridded(self.fn_dat, self.fn_dom, multiple=True, 
                             config=cfg.fn_cfg_nemo)
        
        # Extract latitude and longitude array
        lon = nemo.dataset.longitude.values.squeeze()
        lat = nemo.dataset.latitude.values.squeeze()
        
        # Extract time indices between start and end dates
        t_ind = nemo.dataset.time.values>=self.start_date
        nemo.dataset = nemo.dataset.isel(t_dim=t_ind)
        t_ind = nemo.dataset.time.values<self.end_date
        nemo.dataset = nemo.dataset.isel(t_dim=t_ind)
        
        nemo.dataset = landmask.add_landmask(nemo.dataset)
        nemo.dataset.rename({"depth_0": "depth"})
        
        # Extract model variables.
        nemo.dataset = nemo.dataset[["temperature",
                                          "salinity",
                                          "bathymetry",
                                          "bottom_level",
                                          "landmask"]]
        return nemo.dataset

    def bin_data(self, ds):

        lons = np.linspace(bdy.lonbounds[0], bdy.lonbounds[1],21)
        lats = np.linspace(bdy.latbounds[0], bdy.latbounds[1],21)
        lon_labs = (lons[1:] + lons[:-1]) / 2
        lat_labs = (lats[1:] + lats[:-1]) / 2
        ds = ds.groupby_bins('longitude', lons, labels=lon_labs).mean()
        ds = ds.groupby_bins('latitude', lats, labels=lat_labs).mean()

        return ds

    def grid_en4_profiles(self):
        """
        Create gridded surface en4 data.

        Load en4 data. Expand dimentions to time, lat and lon. Then bin
        over coarse lat/lon intervals.
        """
    
        # create en4 profile object 
        prof = coast.Profile(config=cfg.fn_cfg_prof)
        
        en4_monthly = []
        for date in self.month_list:
            # load en4
            date_str = date.astype("str").replace("-","")
            print (date_str)
            fn_append = "_processed_" + date_str + ".nc"
            en4 = xr.load_dataset(self.en4_path + fn_append)

            # select variables
            en4 = en4[["potential_temperature", "practical_salinity", "depth"]]

            en4 = en4.set_coords("depth")
            en4 = en4.where(en4.depth < 5, drop=True).mean("z_dim")
            en4 = en4.dropna(dim="id_dim", how="all")

            # expand to 3D 
            en4 = en4.set_index(id_dim=["time","longitude","latitude"])
            en4 = en4.unstack("id_dim")

            # bin over lat/lons
            en4 = self.bin_data(en4)

            # time average for each lat/lon bin
            en4 = en4.mean("time").expand_dims("t_dim")

            # assign month coordinate
            time_arr = xr.DataArray([date], dims="t_dim")
            en4 = en4.assign_coords({"time":time_arr})
            en4_monthly.append(en4)

        # join months and save
        en4_all = xr.concat(en4_monthly, dim="t_dim")
        fout_append = "gridded_en4_monthly_mean_%s_%s.nc"%(
                      self.month_list[0], self.month_list[-1])
        en4_all.to_netcdf(cfg.dout_en4 + fout_append)
            
    def get_gridded_en4_climatology(self):
        """ Get seasonal climatololy from gridded en4 data """

        # read gridded surface en4
        fout_append = "gridded_en4_monthly_mean_%s_%s.nc"%(
                      self.month_list[0], self.month_list[-1])
        en4 = xr.open_dataset(cfg.dout_en4 + fout_append)
        
        # calculate climatology
        clim = coast.Climatology()
        clim_mean = clim.make_climatology(en4, "season", 
                      fn_out=self.fn_out + "en4_gridded_surface_climatology.nc")
            
gr = gridded_en4()
#m_ds = gr.model_data()
#o_ds = gr.en4_profiles()
gr.get_gridded_en4_climatology()
