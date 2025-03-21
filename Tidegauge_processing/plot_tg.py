from config import config
import sys

config = config() # initialise variables in python

# IF USING A DEVELOPMENT BRANCH OF COAST, ADD THE REPOSITORY TO PATH:
sys.path.append(config.coast_repo)

import coast
import xarray as xr
import numpy as np
#import datetime
import pandas as pd
import glob, os
from coast import crps_util as cu
import numpy as np

import time
from validate_ssh_tg_hourly import extract_ssh, analyse_ssh, plot_single_cfg, plot_stats_ssh_hourly_compare_cfgs, plot_taylor_tide
import numpy as np

#constit = ['M2','S2','N2','K1','O1','P1','M4']

#args = sys.argv

#port_id = int(args[1])

## Use validate_ssh_tg_hourly.py
if(0):  # Takes 7mins for one month. 1 yr in 2.5 hrs
  fn_nemo_data = config.fn_nemo_data.replace("YYYYMM*_shelftmb", '%04d%02d'%(year,month)+"*_shelftmb")
  fn_extr_out = config.fn_extr_out.replace(".nc", '_%04d%02d'%(year,month)+".nc")
  print(f"Start extract_ssh")
  print(f"Open file, fn_nemo_data:{fn_nemo_data}")
  extract_ssh(fn_nemo_data, config.fn_nemo_domain, config.fn_nemo_cfg,
		config.fn_obs, fn_extr_out,
                chunks = {'time_counter':100}, dist_omit = 5,
		)


if(0):
  fn_extr_out = config.fn_extr_out.replace(".nc", "_??????.nc")
  fn_analyse_out = config.fn_analyse_out.replace(".nc", '_port%02d'%port_id+".nc")
  print(f"Start analyse_ssh")
  analyse_ssh(fn_extr_out, fn_analyse_out, thresholds = np.arange(-.4, 2, 0.1),
                constit_to_save = ['M2','S2','K2','N2','K1','O1','P1','Q1'],
                semidiurnal_constit = ['M2','S2','K2','N2'],
                diurnal_constit = ['K1','O1','P1','Q1'],
                apply_ntr_filter = True,
		port_id = port_id )


if(0):
  for count, file in enumerate(glob.glob(config.fn_analyse_out.replace(".nc", "_port*.nc"))):
    print(file)
    ds = xr.open_dataset(file)
    if count == 0: ds_all = ds
    else:
      ds_all = xr.concat((ds_all, ds), dim='id_dim')
  ds_all.to_netcdf( config.fn_analyse_out )



# Set output file, unless edited fn_ssh_hourly_stats==fn_extr_out
if(1):  # THIS PRODUCES CORR BUT GETS STUCK WITH HARMONICS PLOTS - NO DATA.
  # mkdir config.dn_out
  #plot_single_cfg( config.fn_analyse_out, config.dn_out, config.run_name, file_type='.png')

  # plot modified Taylor Tide diagram
  ref_run_name = "P0.0"
  plot_taylor_tide( [config.fn_analyse_out.replace(config.run_name, ref_run_name), config.fn_analyse_out], config.dn_out, [ref_run_name, config.run_name], file_type='.png')
  
  # plot the differences against a reference run: "diff = run_name - ref_run_name"
  ref_run_name = "P0.0"
  fn_analyse_out_ref = config.fn_analyse_out.replace(config.run_name, ref_run_name)
  plot_stats_ssh_hourly_compare_cfgs( config.fn_analyse_out, fn_analyse_out_ref, config.dn_out, config.run_name, ref_run_name, file_type='.png')
