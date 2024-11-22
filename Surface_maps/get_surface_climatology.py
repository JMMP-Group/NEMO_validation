### work in progress ###
### replication of  EN4_processing/GEN_MOD_Dave_example_profile_validation.py
### purpose > make gridded version of en4 data?

from PythonEnvCfg.config import config, bounds

cfg = config() # initialise variables in python
bdy = bounds("AMM15")

from _utils import landmask
from dask.diagnostics import ProgressBar
import xarray as xr
import pandas as pd
import os
import numpy as np

class surface_climatology(object):

    def __init__(self):

        # File paths 
        self.fn_out = cfg.dn_out + 'surface_maps/'
        
        # Make out directory
        print(os.popen(f"mkdir -p {cfg.dn_out}").read())
        
    def get_season_bias(self):
        """ merge seasonal profile bias files into one dataset """


        ds_list = []
        for season in ['DJF','MAM','JJA','SON']:
            # get data
            path = cfg.dn_out + 'profiles/' + season + "_PRO_DIFF.nc"
            ds = xr.open_dataset(path, chunks='auto')

            ds_list.append(ds.assign_coords(
                           season=('id_dim',[season]*len(ds.id_dim))))

        self.season_ds = xr.concat(ds_list, dim='id_dim')

    def restrict_to_surface(self, depth_lim=5, save=True):
        """ restrict to top x-metres """

        # set bounds
        bounds = [0,depth_lim]

        # sel is faster than where
        self.season_ds = self.season_ds.swap_dims({"z_dim":"depth"})
        self.season_ds = self.season_ds.sel(depth=bounds, method="nearest")

        # average
        self.season_ds = self.season_ds.mean("depth")

    def save_ds(self, ds, path):
        """ save data to path """

        with ProgressBar():
            ds.to_netcdf(self.fn_out + path)

sc = surface_climatology()
sc.get_season_bias()
sc.restrict_to_surface()
sc.save_ds(sc.season_ds, "near_surface_EN4_bias_by_season.nc")
