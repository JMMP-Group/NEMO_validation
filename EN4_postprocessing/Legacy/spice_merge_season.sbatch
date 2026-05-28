#!/bin/bash -l
#SBATCH --mem=20000
#SBATCH --ntasks=1
#SBATCH --output=LOGS/%A_%a.out
#SBATCH --error=LOGS/%A_%a.err
#SBATCH --time=15
source config.sh
conda activate $CONDA_ENV 

echo "time python merge_season.py $1 $2"  # $1 <-- MOD, $2 <-- season
time python merge_season.py $1 $2 > LOGS/merge_season_$1_$2.log