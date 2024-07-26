#!/bin/bash
#SBATCH --partition=short-serial
#SBATCH --mem=40000
#SBATCH -o LOGS/%A_%a.out
#SBATCH -e LOGS/%A_%a.err
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
module add jaspy
source activate $CONDA_ENV

#export RUN_NAME="GS1p1_tide"
#export RUN_NAME="GS1p2_full"
#export RUN_NAME="FES2014"

export RUN_NAME=$1
source config.sh  # sets paths for config given #RUN_NAME

# pass config string
python ./pre_process_harm.py > LOGS/OUT_$1.log
