#!/bin/bash

## Set the machine to be used. Pick one.
export MACHINE="SPICE"  # resource at MO

## Years to loop over during monthly preprocessing: called in iter_en4_proc.sh, iter_sub_METEST.sh
export STARTYEAR=2013 #1980  # 2004
export ENDYEAR=2013 #1980    # 2014

## Process monthly data. Required in iter_sub_METEST.sh
export MOD="P1.5c"  # Model reference name
export GRID="domain_cfg_sf12.nc"  # contains the grid information for NEMO
# options?: "domain_cfg_MEs_01-003_opt_v1.nc" #  "GEG_SF12.nc"
# options?: "GEG_SF12.nc" CO7_EXACT_CFG_FILE.nc


#Conda environment path
export CONDA_ENV="/home/users/o.lambertbrown/.conda/envs/nemo_valid"  ## OR WHAT IS THIS REALLY? WHAT IS THE PATH?

# location of COAsT repo, if using a particular branch
#export COAST_REPO_OLD="/data/users/fred/SET_UP_CONDA_ARTIFACTORY/COAST_SCIPY"
export COAST_REPO_NEW="/home/users/o.lambertbrown/COAsT"  ## OR WHAT IS THIS REALLY? WHAT IS THE PATH?
export COAST_REPO=${COAST_REPO_NEW}
#export COAST_REPO="/data/users/fred/SET_UP_CONDA_ARTIFACTORY/COAST_NOV_2022_DEVELOP/COAsT"
# COAsT configuration files 
export DN_CFG_NEW="/home/users/o.lambertbrown/COAsT/config/"
export DN_CFG=${DN_CFG_NEW}
export FN_CFG_PROF=$DN_CFG"example_en4_profiles.json"
export FN_CFG_NEMO=$DN_CFG"example_nemo_grid_t.json"


# location of raw EN4 data
export DIN_EN4="/data/users/o.lambertbrown/Datasets/EN4_profiles/2013/"
# location of preprocessed EN4 data
export DOUT_EN4="/data/users/o.lambertbrown/Datasets/EN4_profiles/2013/temp"
export REGION="AMM15"  # prefix for preprocessed EN4 data (chunked into files by region and time)

# directory for NEMO domain_cfg.nc file
export DN_DOM="/data/users/o.lambertbrown/Datasets/models/CO9/u-cu674/"
#fn_dom = "/data/users/fred/ME_DOMAINS/
# directory for NEMO data files
#fn_dat = "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/%s/DAILY/%s%02d*T.nc*"%(exper,startyear,month)
export DN_DAT="/data/users/o.lambertbrown/Datasets/models/CO9/u-cu674/"
export DN_OUT="/data/users/o.lambertbrown/Datasets/models/CO9/u-cu674/analysis/"
