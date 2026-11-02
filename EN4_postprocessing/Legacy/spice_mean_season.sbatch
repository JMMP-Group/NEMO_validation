#!/bin/bash -l
#SBATCH --mem=20000
#SBATCH --ntasks=1
#SBATCH --output LOGS/%A_%a_mean.out
#SBATCH --error LOGS/%A_%a_mean.err
#SBATCH --time=30:00
source config.sh
conda activate $CONDA_ENV 

echo "time python mean_season.py $1 $2"  # $1 <-- MOD, $2 <-- season
# NB was ipython
time python mean_season.py $1 $2 > LOGS/mean_season_$1_$2.log
