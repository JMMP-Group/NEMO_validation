import os
  
#os.chdir('/home/users/jelt/dbyrne/code/model_validation/')

from analyse_ssh_hourly import analyse_ssh_hourly, plot_stats_ssh_hourly_single_cfg
from validate_ssh_tg_hourly import extract_ssh, analyse_ssh, plot_single_cfg
import numpy as np

run_name='P0.0'
mass = 'xu-cb676'

#fn_nemo_data = "/gws/nopw/j04/jmmp_collab/CO9_AMM15/outputs/hourly/{0}/*.nc".format(run_name)
fn_nemo_data = "/gws/nopw/j04/jmmp/MASS/{0}/shelftmb/2014*_shelftmb_grid_T.nc".format(mass)
fn_nemo_domain = "/gws/nopw/j04/jmmp_collab/CO9_AMM15/inputs/domains/CO7_EXACT_CFG_FILE.nc"
fn_nemo_cfg = "/home/users/jelt/GitHub/COAsT/config/example_nemo_grid_t.json"

fn_out = "/gws/nopw/j04/jmmp/CO9_AMM15_validation/{0}/tg_analysis/ssh_hourly_{1}.nc".format(run_name.upper(), run_name)
#fn_out = "/home/users/dbyrne/CO9_AMM15/analysis/ssh_hourly_{0}.nc".format(run_name)
fn_obs = "/gws/nopw/j04/jmmp/CO9_AMM15/obs/tg_amm15.nc"

constit = ['M2','S2','N2','K1','O1','P1','M4']

#analyse_ssh_hourly(fn_nemo_data, fn_nemo_domain, fn_nemo_cfg, fn_obs, fn_out, constit_to_save=constit, chunks = {'time_counter':50})

dn_out = "/home/users/jelt/GitHub/NEMO_validation/tidegauge_analysis/FIGS"

## Use validate_ssh_tg_hourly.py
fn_extr_out = "/gws/nopw/j04/jmmp/CO9_AMM15_validation/{0}/tg_analysis/ssh_hourly_extract_{1}.nc".format(run_name.upper(), run_name)
if(0):  # Takes 7mins for one month. 1 yr in 2.5 hrs
  print(f"Start extract_ssh")
  extract_ssh(fn_nemo_data, fn_nemo_domain, fn_nemo_cfg, fn_obs, fn_extr_out,
                     chunks = {'time_counter':100}, dist_omit = 5)

fn_analyse_out = "/gws/nopw/j04/jmmp/CO9_AMM15_validation/{0}/tg_analysis/ssh_hourly_analyse_{1}.nc".format(run_name.upper(), run_name)

if(0):
  print(f"Start analyse_ssh")
  analyse_ssh(fn_extr_out, fn_analyse_out, thresholds = np.arange(-.4, 2, 0.1),
                constit_to_save = ['M2','S2','K2','N2','K1','O1','P1','Q1'],
                semidiurnal_constit = ['M2','S2','K2','N2'],
                diurnal_constit = ['K1','O1','P1','Q1'],
                apply_ntr_filter = True )
if(1):  # THIS PRODUCES CORR BUT GETS STUCK WITH HARMONICS PLOTS - NO DATA.
  fn_ssh_hourly_stats = fn_analyse_out
  #fn_ssh_hourly_stats = "/gws/nopw/j04/jmmp/CO9_AMM15_validation/P0.0/tg_analysis/ssh_hourly_analyse_2014_P0.0.nc"
  plot_single_cfg( fn_ssh_hourly_stats, dn_out, run_name, file_type='.png')

## Use analyse_ssh_hourly.py
fn_analyse_out_b = "/gws/nopw/j04/jmmp/CO9_AMM15_validation/{0}/tg_analysis/ssh_hourly_analyse_b_{1}.nc".format(run_name.upper(), run_name)
if(0):
  analyse_ssh_hourly( fn_nemo_data, fn_nemo_domain, fn_nemo_cfg, fn_obs, fn_analyse_out_b,
             thresholds=np.arange(0, 2, 0.1),
             constit_to_save=['M2', 'S2', 'K1', 'O1'],
             chunks={'time_counter': 100})
if(0):  # THIS PRODUCES PLOTS WITH NaNs? FOR VALUES ON MAPS
  fn_ssh_hourly_stats = fn_analyse_out_b
  plot_stats_ssh_hourly_single_cfg( fn_ssh_hourly_stats, dn_out, run_name, file_type='.png')
