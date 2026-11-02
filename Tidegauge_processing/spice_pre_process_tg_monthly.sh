#!/bin/bash -l
#SBATCH --mem=40000
#SBATCH -o LOGS/%A_%a.out
#SBATCH -e LOGS/%A_%a.err
#SBATCH --time=60:00
#SBATCH --ntasks=1
conda activate $CONDA_ENV

# pass integer year and month
python ./pre_process_tg_monthly.py $1 $2 > LOGS/OUT_$1_$2.log
