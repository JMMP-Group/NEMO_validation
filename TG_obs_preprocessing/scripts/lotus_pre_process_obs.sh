#!/bin/bash
#SBATCH --partition=short-serial
#SBATCH --mem=40000
#SBATCH -o LOGS/%A_%a.out
#SBATCH -e LOGS/%A_%a.err
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
module add jaspy

source config.sh  # sets $MACHINE. Is called again before other variables are used

rm LOGS/OUT* LOGS/*.err LOGS/*.out

source activate $CONDA_ENV

python ./pre_process_obs.py > LOGS/OUT_pre_process_obs.log
