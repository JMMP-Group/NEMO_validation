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
#from validate_ssh_tg_hourly import extract_ssh, analyse_ssh, plot_single_cfg
import numpy as np

#constit = ['M2','S2','N2','K1','O1','P1','M4']

#args = sys.argv

#year = int(args[1])
#month  = int(args[2])

def amp_pha_from_re_im(creal,cimag):
    """
    Example usage: amp,pha = amp_pha_from_re_im(ds.M2x,ds.M2y)
    """
    cc=creal+cimag*1j
    amp=np.abs(cc)
    pha=np.angle(cc)*180/np.pi
    return(amp,pha)

def load_and_save_nemo():
    nemo = coast.Gridded(config.fn_nemo_data, config.fn_nemo_domain, config=config.fn_nemo_cfg, multiple=True)  # , chunks=chunks)

    # nemo.dataset['M2y'] = -nemo.dataset.M2y  # Not particularly happy about this fix...
    # print("nemo.dataset['M2y'] = -nemo.dataset.M2y  # Not particularly happy about this fix... ")
    nemo.dataset['A'], nemo.dataset['G'] = amp_pha_from_re_im(nemo.dataset.M2x, nemo.dataset.M2y)

    # Create a landmask array in Gridded
    nemo.dataset["landmask"] = nemo.dataset.bottom_level == 0
    nemo.dataset = nemo.dataset.rename({"depth_0": "depth"})
    # Find nearest NEMO geographic neighbours to observation locations
    print(f"Implementing the obs_operator: obs.obs_operator(nemo)")
    tg = obs.obs_operator(nemo)
    print(f"Write processed file to {config.fn_analysis_out}")
    tg.to_netcdf(config.fn_analysis_out)
    return tg

def load_and_save_fes2014():
    ## Load FES model harmonics data

    # Load FES data
    fes     = coast.Gridded(config.fn_fes_amp, config.fn_nemo_domain, config=config.fn_nemo_cfg)
    fes_pha = coast.Gridded(config.fn_fes_pha, config.fn_nemo_domain, config=config.fn_nemo_cfg)
    fes.dataset['A'] = fes.dataset.M2amp
    fes.dataset['G'] = fes_pha.dataset.M2pha

    # Create a landmask array in Gridded
    fes.dataset["landmask"] = fes.dataset.bottom_level == 0
    fes.dataset = fes.dataset.rename({"depth_0": "depth"})

    # Find nearest NEMO geographic neighbours to observation locations
    print(f"Implementing the obs_operator: obs.obs_operator(nemo)")
    tg = obs.obs_operator(fes)
    print(f"Write processed file to {config.fn_analysis_out}")
    tg.to_netcdf(config.fn_analysis_out)
    return tg



# Load obs data
###############
obs = coast.Tidegauge(dataset=xr.open_dataset(config.fn_harm_obs))
obs.dataset = obs.dataset.rename_dims({"locs":"id_dim"})  # coast tidegauge object expects dims {id_dim, t_dim}
obs.dataset = obs.dataset.rename_vars({"z1":"M2x", "z2":"M2y"})

#obs.dataset['A'] = np.sqrt(np.square(obs.dataset.M2x) + np.square(obs.dataset.M2y))
#obs.dataset['G'] = -np.arctan2(obs.dataset.M2x, obs.dataset.M2y)  # z1=M2x=a.sin(g); z2=M2y=a.cos(g). Phase increasing _clockwise_ from 'noon'
obs.dataset['A'], obs.dataset['G'] = amp_pha_from_re_im(obs.dataset.M2x, obs.dataset.M2y)

# Load model data
#################
if config.run_name == "fes2014":
    print(f"run_name: {config.run_name}")
    tg = load_and_save_fes2014()
else:
    print(f"run_name: {config.run_name}")
    tg = load_and_save_nemo()

