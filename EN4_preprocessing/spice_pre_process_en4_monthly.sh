#!/bin/bash -l
#SBATCH --mem=200000
#SBATCH --ntasks=1
#SBATCH --output=LOGS/%A_%a.out
#SBATCH --error=LOGS/%A_%a.err
#SBATCH --time=180
source config.sh
conda activate $CONDA_ENV

#python pre_process_en4_1990_2020.py
# Pass in integer year and month: python pre_process_en4_monthly.py <year> <month>
python pre_process_en4_monthly.py $1 $2 > LOGS/OUT_$1_$2.log
