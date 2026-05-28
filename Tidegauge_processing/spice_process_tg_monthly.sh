#!/bin/bash -l
#SBATCH --mem=40000
#SBATCH -o LOGS/%A_%a.out
#SBATCH -e LOGS/%A_%a.err
#SBATCH --time=60:00
#SBATCH --ntasks=1
conda activate $CONDA_ENV

echo $1  # port number
# pass integer port number
python ./process_tg_monthly.py $1 > LOGS/OUT_$1.log
