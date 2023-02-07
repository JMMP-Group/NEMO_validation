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
import os
from coast import crps_util as cu
import numpy as np

import time
from validate_ssh_tg_hourly import extract_ssh, analyse_ssh, plot_single_cfg
import numpy as np

#constit = ['M2','S2','N2','K1','O1','P1','M4']


## Use validate_ssh_tg_hourly.py
if(0):  # Takes 7mins for one month. 1 yr in 2.5 hrs
  print(f"Start extract_ssh")
  extract_ssh(config.fn_nemo_data, config.fn_nemo_domain, config.fn_nemo_cfg,
		config.fn_obs, config.fn_extr_out,
                chunks = {'time_counter':100}, dist_omit = 5)


if(0):
  print(f"Start analyse_ssh")
  analyse_ssh(config.fn_extr_out, config.fn_analyse_out, thresholds = np.arange(-.4, 2, 0.1),
                constit_to_save = ['M2','S2','K2','N2','K1','O1','P1','Q1'],
                semidiurnal_constit = ['M2','S2','K2','N2'],
                diurnal_constit = ['K1','O1','P1','Q1'],
                apply_ntr_filter = True )

if(1):  # THIS PRODUCES CORR BUT GETS STUCK WITH HARMONICS PLOTS - NO DATA.
  plot_single_cfg( config.fn_ssh_hourly_stats, config.dn_out, config.run_name, file_type='.png')
