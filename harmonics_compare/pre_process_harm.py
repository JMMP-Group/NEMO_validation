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

constit_list = ['M2','S2','N2','K1','O1','Q1','P1','K2'] #,'M4','MS4']

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
    for count, constit in enumerate(constit_list):
        # Load FES data
        #fes     = coast.Gridded(config.fn_fes_amp, config.fn_nemo_domain, config=config.fn_nemo_cfg)
        #fes_pha = coast.Gridded(config.fn_fes_pha, config.fn_nemo_domain, config=config.fn_nemo_cfg)
        # Rename coords + dims. Convert to 2d lat/lon for COAsT interpolation to work
        try:
            ds = xr.open_dataset(config.dn_fes+constit+"_z.nc")
            ds = ds.rename_dims({'latitude':'y_dim', 'longitude':'x_dim'})
            #ds_latitude, ds_longitude = np.meshgrid(ds.latitude.values, ds.longitude.values)
            lat, lon = ds.latitude, ds.longitude
            ds_latitude, ds_longitude = xr.broadcast(lat, lon)
            ds = ds.drop_vars(["latitude", "longitude"])
            ds['latitude'] = xr.DataArray(ds_latitude, dims=["y_dim", "x_dim"], coords=[lat, lon])
            ds['longitude'] = xr.DataArray(ds_longitude, dims=["y_dim", "x_dim"], coords=[lat, lon])
            ds = ds.set_coords(['latitude', 'longitude'])
            ds = ds.drop_vars(["x_dim", "y_dim"])
            #print(ds.latitude)
            #ds['latitude'], ds['longitude'] = np.meshgrid(ds.latitude.values, ds.longitude.values)

            if (count == 0) and (constit == "M2"):
              fes = coast.Gridded()
              fes.dataset = ds
              fes.dataset = fes.dataset.rename({"amplitude":"M2amp", "phase":"M2pha"})

            #print(fes.dataset)
            #print(ds)

            #nemo = coast.Gridded(config.fn_nemo_data, config.fn_nemo_domain, config=config.fn_nemo_cfg, multiple=True)  # , chunks=chunks)
            #print(f'\n')
            #print(nemo.dataset)

            #fes.dataset = fes.dataset.rename({"amplitude": "A", "phase": "G"})

            fes.dataset[constit+'x'] = xr.zeros_like(ds.amplitude)
            fes.dataset[constit+'y'] = xr.zeros_like(ds.amplitude)
            fes.dataset[constit+'x'], fes.dataset[constit+'y'] = re_im_from_amp_pha(ds['amplitude'], ds['phase'])

        except:
            print(f"load and save FES2014: Skipped constituent: {constit}")

    # Create a landmask array in Gridded
    #fes.dataset["landmask"] = fes.dataset.bottom_level == 0
    #fes.dataset = fes.dataset.rename({"depth_0": "depth"})

    # Find nearest NEMO geographic neighbours to observation locations
    print(f"Implementing the obs_operator: obs.obs_operator(nemo)")
    tg = obs.obs_operator(fes)
    print(f"Write processed file to {config.fn_analysis_out}")
    tg.dataset.attrs['title'] = "FES2014 complex harmonics at observation locations. (Using https://github.com/JMMP-Group/NEMO_validation)"    
    tg.dataset.attrs['history'] = "Created "+str(np.datetime64('now'))
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
    obs.dataset.to_netcdf(config.fn_analysis_out.replace(config.run_name,'obs'))
else:
    print(f"run_name: {config.run_name}")
    tg = load_and_save_nemo()

