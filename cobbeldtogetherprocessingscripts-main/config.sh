#!/bin/bash

## Years to loop over during monthly preprocessing: called in iter_en4_proc.sh
export STARTYEAR=1980
export ENDYEAR=1980

## Process monthly data
export MOD="P0.0"  # Model name
export GRID="CO7_EXACT_CFG_FILE.nc"  # contains the grid information for NEMO
# options?: "domain_cfg_MEs_01-003_opt_v1.nc" #  "GEG_SF12.nc"
# options?: "GEG_SF12.nc" CO7_EXACT_CFG_FILE.nc

## Set the machine to be used. Pick one.
export MACHINE="LOTUS"  # resource on JASMIN
#export MACHINE="SPICE"  # resource at MO

## SETTINGS FOR JASMIN LOTUS PROCESSING
if [ $MACHINE = "LOTUS" ]; then
  # Use to activate conda environment in $MACHINE_pre_process_en4_monthly.sh
  export CONDA_ENV="/home/users/jelt/miniconda3/envs/coast_dev"

  # location of COAsT repo, if using a particular branch
  export COAST_REPO="/home/users/jelt/GitHub/COAsT"
  # COAsT configuration file for EN4 profiles
  export FN_CFG_PROF="/home/users/jelt/GitHub/COAsT/config/example_en4_profiles.json"

  # location of raw EN4 data
  export DIN_EN4="/home/users/jelt/EN4/"

  # prefix for preprocessed EN4 data (chunked into files by region and time)
  export REGION="AMM15"
  # location of preprocessed EN4 data
  #export DOUT_EN4="/home/users/jelt/tmp/"
  export DOUT_EN4="/home/users/jelt/"$REGION"/"

## SETTINGS FOR MET OFFICE SPICE PROCESSING
elif [ $MACHINE = "SPICE" ]; then
  # Use to activate conda environment in $MACHINE_pre_process_en4_monthly.sh
  export CONDA_ENV_OLD="~/envs/coast"  ## WHAT IS THIS REALLY? THIS IS YOUR OLD COAST ENV
  # Use to active conda environment in the "process monthly data" model/en4 inter-comparison
  export CONDA_ENV_NEW="coast_nov2022"  ## OR WHAT IS THIS REALLY? WHAT IS THE PATH?
  export CONDA_ENV=CONDA_ENV_OLD  ## Ideally you would only ever use one (new) conda environment

  # location of COAsT repo, if using a particular branch
  export COAST_REPO = "/data/users/fred/SET_UP_CONDA_ARTIFACTORY/COAST_SCIPY"
  # COAsT configuration file for EN4 profiles
  export FN_CFG_PROF="/data/users/fred/coast_demo/config/example_en4_profiles.json"

  # location of raw EN4 data
  export DIN_EN4="/scratch/fred/EN4/"
  # location of preprocessed EN4 data
  export DOUT_EN4="/scratch/fred/EN4/"
  export REGION="SCIPY"  # prefix for preprocessed EN4 data (chunked into files by region and time)


fi
