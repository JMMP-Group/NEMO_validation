#!/usr/bin/env python3
'''
For processing/writing obs for model assessment.
Uses process_obs routine (to be checked for details)

Multi stage process requiring NEMO domain configuration, FES2014 data dir

env: gtm_obs_env

micromamba create -n gtm_obs_env python=3.10 numpy cmocean netCDF4 scipy cartopy matplotlib
micromamba activate gtm_obs_env

jelt 2025-06-04

Initial version:
David Byrne, Jeff Polton & Colin Bell (2021): Creation of a global tide analysis
dataset: Application of NEMO and an offline objective analysis scheme, Journal of Operational
Oceanography, https://doi.org/10.1080/1755876X.2021.2000249
Code repository: https://github.com/NOC-MSM/global_tide


Changelog:
2020-XX-XX: Initial version  Byrne et al. (2021) 
2025-06-04: Ported to NEMO_validation, to be used for model assessment, with minor edits. (jelt)
'''

from config import config
import sys

config = config() # initialise variables in python

# IF USING A DEVELOPMENT BRANCH OF COAST, ADD THE REPOSITORY TO PATH:
sys.path.append(config.coast_repo)


import numpy as np
import GTM_assimilation as ass
import dbp_read as dbr
import dbp_gtm as gtm
import dbp_general as dbg

#const=['M2']
#const = ['SA','SSA','MF','MM','MSF','Q1','O1','P1','S1','K1','J1','2N2','MU2',
#         'N2','NU2','M2','MKS2','L2','T2','S2','R2','K2','M3','MN4','M4','MS4',
#         'S4','M6','M8']  
const = ['Q1','O1','P1','S1','K1','J1',
         'N2','M2','S2','K2',
        ]

# File directory for ensemble arrays, nemo arrays.
if(0): # GTM processing
    fn_nemo_harm      = '/projectsa/NOCglobaltide/NEMO/output/GO8_R12_v1.0_final.nc'
    fn_nemo_dom       = '/projectsa/NOCglobaltide/NEMO/INPUTS/domains/domcfg_eORCA12_10levels_23m.nc'
    file_fes_dir = '/Users/jelt/DATA/FES2014/ocean_tide_extrapolated/'
    #fn_out_dir = '/projectsa/NOCglobaltide/data/obs/for_DA_sparse/obs_'
    fn_out_dir = '/scratch/jelt/GTM_tmp/data/obs/for_DA_sparse/obs_'


fn_nemo_harm      = '/Users/jelt/Downloads/SENEMO/TIDE/SENEMO_1y_19810101_19811231_grid_T_2D.nc'
fn_nemo_dom       = '/Users/jelt/Downloads/SENEMO/TIDE/domain_cfg.nc'
file_fes_dir = '/Users/jelt/DATA/FES2014/ocean_tide_extrapolated/'
fn_out_dir = '/Users/jelt/Downloads/SENEMO/data/obs/for_validation_sparse/obs_'

nemo_stride = 1

n_const = len(const)
for ii in range(0,n_const):
    # Define temporary constituent variable.
    cc = const[ii]
    print(cc)
    tmp = ass.read_obs([cc], dir='/Users/jelt/GitHub/NEMO_validation/TG_obs_preprocessing/data/obs/')
    y_lon = tmp[0]; y_lat = tmp[1]; y_z1 = tmp[2]; y_z2 = tmp[3];
    y_a = tmp[4]; y_g = tmp[5]; obs_id = tmp[6] 
    
    print(obs_id[100])

    # Read FES
    tmp = dbr.read_fes_2D_harm(file_fes_dir, [cc], istride=2, jstride=2)
    fes_lon = tmp[2]; fes_lat = tmp[3]; 
    fes_mask = tmp[4]; fes_z1 = tmp[5]; fes_z2 = tmp[6]
    fes_z1 = np.squeeze(fes_z1); fes_z2 = np.squeeze(fes_z2)
    
    del tmp
    if ii == 0:
        # Read NEMO grid and mask data
        tmp = dbr.read_nemo_2D_harm(fn_nemo_harm, [cc], fn_nemo_dom,
                                        istride = nemo_stride, jstride= nemo_stride,
                                        var='z')
        nemo_lon = tmp[2]; nemo_lat = tmp[3]; nemo_mask = tmp[4]; 
        y_ii, y_jj = dbg.geo_nn_indices(nemo_lon, nemo_lat, y_lon, y_lat, 
                                        mask = nemo_mask)
        fes_ii, fes_jj = dbg.geo_nn_indices(fes_lon, fes_lat, y_lon, y_lat, 
                                    mask = fes_mask)
    
    y_lonP1, y_latP1, yP1, yid1 = ass.process_obs(y_lon, y_lat, y_z1, obs_id, 
                                         nemo_lon, nemo_lat, nemo_mask, 
                                         fes_lon, fes_lat, fes_z1, fes_mask, 
                                         sobs_rad = 200, y_ii = y_ii, y_jj = y_jj,
                                         fes_ii = fes_ii, fes_jj = fes_jj)
    y_lonP2, y_latP2, yP2, yid2 = ass.process_obs(y_lon, y_lat, y_z2, obs_id, 
                                         nemo_lon, nemo_lat, nemo_mask, 
                                         fes_lon, fes_lat, fes_z2, fes_mask, 
                                         sobs_rad = 200, y_ii = y_ii, y_jj = y_jj,
                                         fes_ii = fes_ii, fes_jj = fes_jj)
    
    fn_out = fn_out_dir + cc + '.nc'
    gtm.GTM_file_create_obs(fn_out, len(y_lonP1), cc)
    gtm.GTM_file_append_obs(fn_out, y_lonP1, 'longitude', 'deg')
    gtm.GTM_file_append_obs(fn_out, y_latP1, 'latitude', 'deg')
    gtm.GTM_file_append_obs(fn_out, yP1, 'z1')
    gtm.GTM_file_append_obs(fn_out, yP2, 'z2')
    gtm.GTM_file_append_obs(fn_out, yid1, 'obs_id1')
    gtm.GTM_file_append_obs(fn_out, yid2, 'obs_id2')
