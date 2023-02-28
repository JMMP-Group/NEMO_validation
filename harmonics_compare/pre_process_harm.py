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

def amp_pha_from_re_im(creal, cimag):
    """
    Example usage: amp,pha = amp_pha_from_re_im(ds.M2x, ds.M2y)
    """
    cc=creal+cimag*1j
    amp=np.abs(cc)
    pha=np.angle(cc)*180/np.pi
    return amp, pha

def re_im_from_amp_pha(amp, pha):
    """
    Assumes phase in degrees
    Example usage: amp,pha = amp_pha_from_re_im(ds.amp, ds.pha)
    """
    if np.max(np.abs(pha)) <= 2*np.pi:
        print(f"Warning. Check phase units. Expected degrees. Max: {np.max(pha)}. Min: {np.min(pha)}")
    re = amp*np.cos(pha*np.pi/180)
    im = amp*np.sin(pha*np.pi/180)
    return re, im

def load_and_save_nemo():
    nemo = coast.Gridded(config.fn_nemo_data, config.fn_nemo_domain, config=config.fn_nemo_cfg, multiple=True)  # , chunks=chunks)

    # nemo.dataset['M2y'] = -nemo.dataset.M2y  # Not particularly happy about this fix...
    # print("nemo.dataset['M2y'] = -nemo.dataset.M2y  # Not particularly happy about this fix... ")
    nemo.dataset['A'] = xr.zeros_like(nemo.dataset.M2x)
    nemo.dataset['G'] = xr.zeros_like(nemo.dataset.M2x)
    nemo.dataset['A'].values, nemo.dataset['G'].values = amp_pha_from_re_im(nemo.dataset.M2x, nemo.dataset.M2y)

    # Create a landmask array in Gridded
    nemo.dataset["landmask"] = nemo.dataset.bottom_level == 0
    nemo.dataset = nemo.dataset.rename({"depth_0": "depth"})
    # Find nearest NEMO geographic neighbours to observation locations
    print(f"Implementing the obs_operator: obs.obs_operator(nemo)")
    tg = obs.obs_operator(nemo)
    print(f"Write processed file to {config.fn_analysis_out}")
    tg.dataset.to_netcdf(config.fn_analysis_out)
    return tg

def load_and_save_fes2014():
    ## Load FES model harmonics data

    # Load FES data
    fes     = coast.Gridded(config.fn_fes_amp, config.fn_nemo_domain, config=config.fn_nemo_cfg)
    fes_pha = coast.Gridded(config.fn_fes_pha, config.fn_nemo_domain, config=config.fn_nemo_cfg)
    fes.dataset = fes.dataset.rename({"M2amp": "A"})
    fes.dataset['G'] = fes_pha.dataset.M2pha

    fes.dataset['M2x'] = xr.zeros_like(fes.dataset.A)
    fes.dataset['M2y'] = xr.zeros_like(fes.dataset.A)
    fes.dataset['M2x'], fes.dataset['M2y'] = re_im_from_amp_pha(fes.dataset['A'], fes.dataset['G'])

    # Create a landmask array in Gridded
    fes.dataset["landmask"] = fes.dataset.bottom_level == 0
    fes.dataset = fes.dataset.rename({"depth_0": "depth"})

    # Find nearest NEMO geographic neighbours to observation locations
    print(f"Implementing the obs_operator: obs.obs_operator(nemo)")
    tg = obs.obs_operator(fes)
    print(f"Write processed file to {config.fn_analysis_out}")
    tg.dataset.to_netcdf(config.fn_analysis_out)
    return tg



# Load obs data
###############
obs = coast.Tidegauge(dataset=xr.open_dataset(config.fn_harm_obs))
obs.dataset = obs.dataset.rename_dims({"locs":"id_dim"})  # coast tidegauge object expects dims {id_dim, t_dim}
obs.dataset = obs.dataset.rename_vars({"z1":"M2x", "z2":"M2y"})

#obs.dataset['A'] = np.sqrt(np.square(obs.dataset.M2x) + np.square(obs.dataset.M2y))
#obs.dataset['G'] = -np.arctan2(obs.dataset.M2x, obs.dataset.M2y)  # z1=M2x=a.sin(g); z2=M2y=a.cos(g). Phase increasing _clockwise_ from 'noon'
obs.dataset['A'] = xr.zeros_like(obs.dataset.M2x)
obs.dataset['G'] = xr.zeros_like(obs.dataset.M2x)
obs.dataset['A'].values, obs.dataset['G'].values = amp_pha_from_re_im(obs.dataset.M2x, obs.dataset.M2y)

# Load and save model data
##########################
if config.run_name.upper() == "FES2014":
    print(f"run_name: {config.run_name}")
    tg = load_and_save_fes2014()
    # Save obs data only when processing FES2014 data
    obs.dataset.to_netcdf(config.fn_analysis_out.replace(run_name,'obs'))
else:
    print(f"run_name: {config.run_name}")
    tg = load_and_save_nemo()

