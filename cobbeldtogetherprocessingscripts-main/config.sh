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

if [ $MACHINE = "LOTUS" ]; then
  # Use to activate conda environment in $MACHINE_pre_process_en4_monthly.sh
  CONDA_ENV="/home/users/jelt/miniconda3/envs/coast_dev"

elif [ $MACHINE = "SPICE" ]; then
  # Use to activate conda environment in $MACHINE_pre_process_en4_monthly.sh
  CONDA_ENV_OLD="~/envs/coast"  ## WHAT IS THIS REALLY? THIS IS YOUR OLD COAST ENV
  # Use to active conda environment in the "process monthly data" model/en4 inter-comparison
  CONDA_ENV_NEW="coast_nov2022"  ## OR WHAT IS THIS REALLY? WHAT IS THE PATH?
  CONDA_ENV=CONDA_ENV_OLD  ## Ideally you would only ever use one (new) conda environment
fi
