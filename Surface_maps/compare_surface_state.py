from config import config

config = config() # initialise variables in python

import coast
import xarray as xr
import os

class extract_surface(object):
    def __init__(self):

        # paths
        self.fn_dom = config.dn_dom + config.grid_nc
        self.fn_dat = config.dn_out + "profiles/gridded*.nc"
        self.fn_out = config.dn_out + 'surface_maps/'

    def surface_state_climatology(self):
        ds = xr.open_mfdataset(self.fn_dat, combine="nested", 
                               concat_dim="t_dim", parallel=True)
        clim = coast.Climatology()
        clim_mean = clim.make_climatology(ds, "season",
                       fn_out=self.fn_out + 'surface_state_climatology.nc')

es = extract_surface()
es.surface_state_climatology()
