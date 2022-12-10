#!/bin/bash

## Set the machine to be used. Pick one.
export MACHINE="SPICE"  # resource at MO

## Years to loop over during monthly preprocessing: called in iter_en4_proc.sh, iter_sub_METEST.sh
export STARTYEAR=2004 #1980  # 2004
export ENDYEAR=2005 #1980    # 2014

## Process monthly data. Required in iter_sub_METEST.sh
export MOD="P0.0"  # Model reference name
export GRID="CO7_EXACT_CFG_FILE.nc"  # contains the grid information for NEMO
# options?: "domain_cfg_MEs_01-003_opt_v1.nc" #  "GEG_SF12.nc"
# options?: "GEG_SF12.nc" CO7_EXACT_CFG_FILE.nc


# Use to activate conda environment in $MACHINE_pre_process_en4_monthly.sh
export CONDA_ENV_OLD="/home/h01/fred/NOTES/SET_UP_CONDA_ARTIFACTORY/COAST_SCIPY"  ## WHAT IS THIS REALLY? THIS IS YOUR OLD COAST ENV
# Use to active conda environment in the "process monthly data" model/en4 inter-comparison
export CONDA_ENV_NEW="/data/users/fred/SET_UP_CONDA_ARTIFACTORY/COAST_NOV_2022_DEVELOP/COAsT"  ## OR WHAT IS THIS REALLY? WHAT IS THE PATH?
export CONDA_ENV=CONDA_ENV_OLD  ## Ideally you would only ever use one (new) conda environment

# location of COAsT repo, if using a particular branch
export COAST_REPO = "/data/users/fred/SET_UP_CONDA_ARTIFACTORY/COAST_SCIPY"
# COAsT configuration files 
export DN_CFG="/data/users/fred/coast_demo/config/"
export FN_CFG_PROF=$DN_CFG"example_en4_profiles.json"
export FN_CFG_NEMO=$DN_CFG"example_nemo_grid_t.json"


# location of raw EN4 data
export DIN_EN4="/scratch/fred/EN4/"
# location of preprocessed EN4 data
export DOUT_EN4="/scratch/fred/EN4/"
export REGION="SCIPY"  # prefix for preprocessed EN4 data (chunked into files by region and time)

# directory for NEMO domain_cfg.nc file
export DN_DOM="/data/users/fred/ME_DOMAINS/"
#fn_dom = "/data/users/fred/ME_DOMAINS/
# directory for NEMO data files
#fn_dat = "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/DAILY/%s%02d*T.nc*"%(exper,startyear,month)
export DN_DAT="/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/"$MOD"/DAILY/"
export DN_OUT="/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/"$MOD"/analysisb/"
