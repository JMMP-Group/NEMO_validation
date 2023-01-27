import os
  
os.chdir('/home/users/jelt/dbyrne/code/model_validation/') 

from analyse_ssh_hourly import analyse_ssh_hourly
from validate_ssh_tg_hourly import extract_ssh

run_name='P0.0'
mass = 'xu-cb676'

#fn_nemo_data = "/gws/nopw/j04/jmmp_collab/CO9_AMM15/outputs/hourly/{0}/*.nc".format(run_name)
fn_nemo_data = "/gws/nopw/j04/jmmp/MASS/{0}/shelftmb/20140101_shelftmb_grid_T.nc".format(mass)
fn_nemo_domain = "/gws/nopw/j04/jmmp_collab/CO9_AMM15/inputs/domains/CO7_EXACT_CFG_FILE.nc"
fn_nemo_cfg = "/home/users/jelt/GitHub/COAsT/config/example_nemo_grid_t.json"

fn_out = "/gws/nopw/j04/jmmp/CO9_AMM15_validation/{0}/tg_analysis/ssh_hourly_{1}.nc".format(run_name.upper(), run_name)
#fn_out = "/home/users/dbyrne/CO9_AMM15/analysis/ssh_hourly_{0}.nc".format(run_name)
fn_obs = "/gws/nopw/j04/jmmp/CO9_AMM15/obs/tg_amm15.nc"

constit = ['M2','S2','N2','K1','O1','P1','M4']

#analyse_ssh_hourly(fn_nemo_data, fn_nemo_domain, fn_nemo_cfg, fn_obs, fn_out, constit_to_save=constit, chunks = {'time_counter':50})

fn_ext_out = "/gws/nopw/j04/jmmp/CO9_AMM15_validation/{0}/tg_analysis/ssh_hourly_extract_{1}.nc".format(run_name.upper(), run_name)
extract_ssh(fn_nemo_data, fn_nemo_domain, fn_nemo_cfg, fn_ext_obs, fn_out,
                     chunks = {'time_counter':100}, dist_omit = 5)