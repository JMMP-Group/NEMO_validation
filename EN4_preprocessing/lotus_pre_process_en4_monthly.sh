#!/bin/bash 
#SBATCH --partition=short-serial 
#SBATCH -o LOGS/%A_%a.out
#SBATCH -o LOGS/%A_%a.err
#SBATCH --time=30:00
#SBATCH --ntasks=1
module add jaspy
source config.sh
source activate $CONDA_ENV

# Pass in integer year and month: python pre_process_en4_monthly.py <year> <month>
python pre_process_en4_monthly.py $1 $2 > LOGS/OUT_$1_$2.log

