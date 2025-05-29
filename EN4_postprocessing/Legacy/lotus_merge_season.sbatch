#!/bin/bash
#SBATCH --partition=short-serial
#SBATCH --mem=20000
#SBATCH -o LOGS/%A_%a.out
#SBATCH -e LOGS/%A_%a.err
#SBATCH --time=15:00
#SBATCH --ntasks=1
module add jaspy
source activate $CONDA_ENV

echo "python merge_season.py $1 $2"  # $1 <-- MOD, $2 <-- season
python merge_season.py $1 $2 > LOGS/merge_season_$1_$2.log
