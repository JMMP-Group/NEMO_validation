#!/bin/bash

## Years to loop over during monthly preprocessing: called in iter_en4_proc.sh
export STARTYEAR=1980
export ENDYEAR=1980

## Set the machine to be used
export MACHINE="LOTUS"  # resource on JASMIN
#export MACHINE="SPICE"  # resource at MO

if [ $MACHINE = "LOTUS" ]; then
  # Use to activate conda environment in $MACHINE_pre_process_en4_monthly.sh
  CONDA_ENV="/home/users/jelt/miniconda3/envs/coast_dev"

elif [ $MACHINE = "SPICE" ]; then
  # Use to activate conda environment in $MACHINE_pre_process_en4_monthly.sh
  CONDA_ENV="~/envs/coast/bin/activate"
fi
