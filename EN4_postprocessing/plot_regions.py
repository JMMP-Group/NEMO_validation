'''
For plotting the analysis regions
'''

from config import config
config = config() # initialise variables in python

config.dn_out = "/Users/jelt/Downloads/"

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append(config.coast_repo)
import coast

#%% File settings
fn_dom_nemo = "%s%s"%(config.dn_dom, config.grid_nc)
fn_cfg_nemo = config.fn_cfg_nemo

print('Doing regional analysis..')

# load nemo lat/lon grid to define regions (as function of bathymetry)
nemo = coast.Gridded(fn_domain=fn_dom_nemo, config=fn_cfg_nemo)

bath = nemo.dataset.bathymetry.values.squeeze()
lon = nemo.dataset.longitude.values.squeeze()
lat = nemo.dataset.latitude.values.squeeze()


# Define Regional Masks

mm = coast.MaskMaker()
masks_list = []
# Add regional mask for whole domain
masks_list.append(np.ones(lon.shape))  # 0
masks_list.append(mm.region_def_nws_north_north_sea(lon, lat, bath))  # 1
masks_list.append(mm.region_def_nws_outer_shelf(lon, lat, bath))  # 2
masks_list.append(mm.region_def_nws_english_channel(lon, lat, bath))  # 3
masks_list.append(mm.region_def_nws_norwegian_trench(lon, lat, bath))  # 4
masks_list.append(mm.region_def_nws_kattegat(lon, lat, bath))  # 5
masks_list.append(mm.region_def_nws_fsc(lon, lat, bath))  # 6
#masks_list.append(mm.region_def_nws_shelf_break(lon,lat,bath))  # 7
masks_list.append(mm.region_def_nws_south_north_sea(lon, lat, bath))  # 7
masks_list.append(mm.region_def_nws_off_shelf(lon, lat, bath))  # 8
masks_list.append(mm.region_def_nws_irish_sea(lon, lat, bath))  # 9

masks_names = ['whole_domain', 'northern_north_sea','outer_shelf','eng_channel','nor_trench',
	    'kattegat', 'fsc', 'southern_north_sea', 'off_shelf', 'irish_sea' ]
print("Size of names is ", len(masks_names[:]))

mask_xr = mm.make_mask_dataset(lon, lat, masks_list, masks_names)


mm.quick_plot(mask_xr)
plt.contourf(lon,lat,bath, levels=(0,10))
plt.savefig("FIGS/maskmaker_quick_plot.png")
plt.close("all")
