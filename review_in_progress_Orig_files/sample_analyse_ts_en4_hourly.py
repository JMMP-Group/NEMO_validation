'''
This is an example script for doing a comparison of hourly NEMO temperature and salinity 
data against EN4 profile data.

The result of this script will be two analysis files. The first contains profile-by-profile
errors/CRPS at the surface and the bottom. The second contains the analysis in the first
averaged into user defined regional boxes.
'''

import os
import datetime
# ADD path to scripts directory for import
os.chdir('/work/n01/n01/dbyrne/CO9_AMM15/code/model_validation/')

import sys
# If using development version of COAsT uncomment:
sys.path.append('/work/n01/n01/dbyrne/CO9_AMM15/code/COAsT')
import coast

import xarray as xr

# Import validation scripts
import validate_ts_en4_hourly


# Name of the run
run_name='co7'

# File paths
fn_nemo_data = "/work/n01/n01/dbyrne/CO9_AMM15/outputs/co7/concatenated/*.nc"
fn_nemo_domain = "/work/n01/n01/dbyrne/CO9_AMM15/inputs/co7/CO7_EXACT_CFG_FILE.nc"
fn_en4 = "/work/n01/n01/dbyrne/CO9_AMM15/obs/en4/*.nc"
fn_coast = "/work/n01/n01/dbyrne/CO9_AMM15/code/COAsT"
fn_out_analysis = "/work/n01/n01/dbyrne/CO9_AMM15/analysis/ts_hourly_{0}_2.nc".format(run_name)
fn_out_regional = "/work/n01/n01/dbyrne/CO9_AMM15/analysis/ts_hourly_{0}_2.nc".format(run_name)

# Open domain_cfg to get latitude, longitude and bathy information for 
# defining regional masks
dom = xr.open_dataset(fn_nemo_domain, chunks={})

# Use COAsT to make predefined AMM regions.
mm = coast.MASK_MAKER()
regional_masks = []
lon = dom.glamt.values.squeeze()
lat = dom.gphit.values.squeeze()
bath = dom.hbatt.rename({'x':'x_dim','y':'y_dim'}).squeeze().load()
regional_masks.append(mm.region_def_nws_north_sea(lon,lat,bath.values))
regional_masks.append(mm.region_def_nws_outer_shelf(lon,lat,bath.values))
regional_masks.append(mm.region_def_nws_english_channel(lon,lat,bath.values))
regional_masks.append(mm.region_def_nws_norwegian_trench(lon,lat,bath.values))
region_names = ['north_sea','outer_shelf','eng_channel','nor_trench']
dom.close()

# Extract and do profile-by-profile analysis of data
validate_ts_en4_hourly.analyse_hourly_en4(fn_nemo_data, fn_nemo_domain, fn_en4, fn_out, 
					  bathymetry=bath, nemo_chunks = {'time_counter':100})

# Do regional analysis
validate_ts_en4_hourly.analyse_regional(fn_out_analysis, fn_nemo_domain, fn_out_regional,
										regional_masks, region_names)
										
