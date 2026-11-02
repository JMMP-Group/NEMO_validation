#!/bin/bash
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --account=jmmp
#SBATCH --ntasks=1
#SBATCH -o LOGS/%A_%a_mean.out
#SBATCH -e LOGS/%A_%a_mean.err
#SBATCH --time=30:00
module add jaspy
source activate $CONDA_ENV

echo "python merge_mean_surface_crps.py $1"  # $1 <-- MOD
python merge_mean_surface_crps.py $1 > LOGS/merge_mean_surface_crps_$1.log
