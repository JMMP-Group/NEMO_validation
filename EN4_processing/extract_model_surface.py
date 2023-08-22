#!/usr/bin/env python3

from config import config
import sys

config = config() # initialise variables in python

from NEMO_validation._utils import landmask
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
        # get data
        ds = coast.Gridded(fn_data=self.fn_dat, fn_domain=self.fn_dom,
                           multiple=True, config=config.fn_cfg_nemo).dataset

        # rename depth 
        ds = ds.rename({"depth_0": "depth"})

        # select variables
        ds = ds[["temperature","salinity","bathymetry","bottom_level"]]

        # add land mask 
        ds = landmask.add_landmask(ds)

        # mask below surface layer 
        surf_mask = ds.where(ds.depth <= self.surf_lim)
        
        # mean over surface
        extracted = surf_mask.mean(dim="z_dim")

        # save
        extracted = extracted.assign_attrs(dict(case=config.case)) 
        self.fend = "profiles/gridded_model_surface_data_%d%02d.nc"%(self.year,
                                                                     self.month)
        extracted.to_netcdf(config.dn_out + self.fend)

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

