#!/bin/bash
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --mem=40000
#SBATCH -o LOGS/%A_%a.out
#SBATCH -e LOGS/%A_%a.err
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --account=class
module add jaspy
source activate $CONDA_ENV

source config.sh  # sets paths for config

python ./pre_process_obs.py > LOGS/OUT_pre_process_obs.log

