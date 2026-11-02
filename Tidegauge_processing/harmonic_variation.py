"""
WORK IN PROGRESS
Harmonic analysis as a function of window length.
Investigation of how confident should you be of an M2 amplitude taken from a 30 day observational record as an unknown
phase of the nodal cycle.

module add jaspy
export CONDA_ENV="/home/users/jelt/miniconda3/envs/coast_dev"
source activate $CONDA_ENV
"""

# location of COAsT repo, if using a particular branch
coast_repo = "/home/users/jelt/GitHub/COAsT"
coast_repo = "/Users/jelt/GitHub/COAsT"

import sys
# IF USING A DEVELOPMENT BRANCH OF COAST, ADD THE REPOSITORY TO PATH:
sys.path.append(coast_repo)

import coast
from coast import general_utils as gu
from coast import plot_util as pu
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from datetime import datetime, timedelta
import pandas as pd
from dateutil.relativedelta import relativedelta
import matplotlib.dates as mdates
import utide as ut
import scipy.signal as signal
import os
import glob
from dask.diagnostics import ProgressBar



fn_obs = "/gws/nopw/j04/jmmp/CO9_AMM15/obs/tg_amm15.nc"
fn_obs = "/Users/jelt/Downloads/tg_amm15.nc"

constit = ['M2','S2','N2','K1','O1','P1','M4']

#analyse_ssh_hourly(fn_nemo_data, fn_nemo_domain, fn_nemo_cfg, fn_obs, fn_out, constit_to_save=constit, chunks = {'time_counter':50})

dn_out = "/home/users/jelt/GitHub/NEMO_validation/tidegauge_analysis/FIGS"
dn_out = "/Users/jelt/GitHub/NEMO_validation/tidegauge_analysis/FIGS" ## JEFF

min_datapoints = 8760  # 12 month of hrly data
constit_to_save = ['M2','S2','K2','N2','K1','O1','P1','Q1']

# Create the object and then inset the netcdf dataset
ds_ssh = coast.Tidegauge(dataset=xr.open_dataset(fn_obs))
try: ds_ssh.dataset = ds_ssh.dataset.swap_dims({"port": "id_dim"})
except: pass
try: ds_ssh.dataset = ds_ssh.dataset.swap_dims({"time": "t_dim"})
except: pass
try: ds_ssh.dataset = ds_ssh.dataset.drop_vars("bad_flag")
except: pass

# Define Dimension Sizes
n_port = ds_ssh.dataset.dims['id_dim']
n_time = 8760 # ds_ssh.dataset.dims['t_dim']
analysis_length = np.arange(25,n_time,1000)  # indices along t_dim dimension. For hourly data this represents hours also
n_length = len(analysis_length)
n_constit = len(constit_to_save)

# ANALYSIS dataset
ds_stats = xr.Dataset(coords=dict(
    longitude=('port', ds_ssh.dataset.longitude.values),
    latitude=('port', ds_ssh.dataset.latitude.values),
    constituent=('constituent', constit_to_save),
    length=('length', analysis_length)),
    data_vars=dict(
        a_obs=(['port', 'constituent', 'length'], np.zeros((n_port, n_constit, n_length)) * np.nan)))

print(ds_ssh.dataset)
print(ds_ssh.dataset.mean(dim="t_dim"))
# Commence analysis
###################
tganalysis = coast.TidegaugeAnalysis()
ds = tganalysis.demean_timeseries(ds_ssh.dataset)

for count, t_dim_ind in enumerate(analysis_length):
    # Harmonic analysis
    ha_obs = tganalysis.harmonic_analysis_utide(ds.dataset.ssh.isel(t_dim=slice(0,t_dim_ind)), min_datapoints=t_dim_ind)

    # save constituents (name + amp/pha)
    for pp in range(n_port):
        if ha_obs[pp] == []: # empty analysis for port
            continue
        else:
            a_dict_obs = dict(zip(ha_obs[pp].name, ha_obs[pp].A))
            for cc in range(0, len(constit_to_save)):
                if constit_to_save[cc] in ha_obs[pp].name:
                    ds_stats['a_obs'][pp, cc, count] = a_dict_obs[constit_to_save[cc]]

print(ds_stats)

# Look at M2 analysis
ds_stats.a_obs.isel(constituent=0).plot()
plt.show()

# A look at the port locations
plt.plot( (ds_ssh.dataset.longitude+180)% 360 - 180, ds_ssh.dataset.latitude, '.')
for i in np.arange(n_port):
    try:
        plt.text(float((ds_ssh.dataset.longitude[i].values+180)% 360 - 180), float(ds_ssh.dataset.latitude[i].values), str(i))
    except:
        print(float((ds_ssh.dataset.longitude[i].values+180)% 360 - 180), float(ds_ssh.dataset.latitude[i].values), str(i))
plt.show()