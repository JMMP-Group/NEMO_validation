#!/usr/bin/env python3

from config import config
import sys

config = config() # initialise variables in python

# IF USING A DEVELOPMENT BRANCH OF COAST, ADD THE REPOSITORY TO PATH:
sys.path.append(config.coast_repo)

import coast
import xarray as xr
import numpy as np
import os
from coast import crps_util as cu
import time

class extract_surface(object):
    """ 
    Extract surface data from model

    To run: e.g.
    python surface_crps.py 2004 1
    """

    def __init__(self, year, month, surf=5):

        # surface averaging depth
        self.surf_lim = surf 

        # assign dates
        self.year = year
        self.month = month

        # paths
        self.fn_dom = config.dn_dom + config.grid_nc
        self.fn_dat = "%s%s%02d*T.nc"%(config.dn_dat, year, month)

        # ensure directories exist
        print(os.popen(f"mkdir -p {config.dn_out}").read())

    def extract(self):
        """
        extract NEMO surface data 
        """
        # create nemo object
        nemo = coast.Gridded(fn_data=self.fn_dat, fn_domain=self.fn_dom,
                             multiple=True, config=config.fn_cfg_nemo)

        # assign case number
        nemo.dataset = nemo.dataset.assign_attrs(dict(case=config.case)) 
        
        # rename depth 
        nemo.dataset = nemo.dataset.rename({"depth_0": "depth"})

        # select variables
        nemo.dataset = nemo.dataset[["temperature", 
                                     "salinity",
                                     "bathymetry",
                                     "bottom_level"]]

        # add land mask
        nemo.dataset = self.add_landmask(nemo.dataset)

        # mask below surface layer 
        surf_mask = nemo.dataset.where(nemo.dataset.depth <= self.surf_lim)
        
        # mean over surface
        self.extracted = surf_mask.mean(dim="z_dim")

        # save
        self.fend = "gridded_model_surface_data_%d%02d.nc"%(self.year,
                                                            self.month)
        self.extracted.to_netcdf(config.dn_out + self.fend)

    def add_landmask(self, ds): 
        """
        Create a landmask array -- important for obs_operator. Calculated 
        from bottom_level.
        """

        # add landmask
        ds["landmask"] = ds.bottom_level == 0

        return ds

if __name__ == '__main__':
    #import dask
    #dask.config.set(scheduler='synchronous') 

    # user input
    args = sys.argv
    year = int(args[1])
    month = int(args[2])

    extr = extract_surface(year, month)
    extr.extract()

    print(f"""Surface data extracted and written to file:
              {config.dn_out+extr.fend}
              Ready for surface or crps processing
              e.g. merge_mean_surface_crps.py""")

