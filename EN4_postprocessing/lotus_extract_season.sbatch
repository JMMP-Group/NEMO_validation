#!/bin/bash
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --account=jmmp
#SBATCH --mem=20000
#SBATCH -o LOGS/%A_%a.out
#SBATCH -e LOGS/%A_%a.err
#SBATCH --time=15:00
#SBATCH --ntasks=1
module add jaspy
source activate $CONDA_ENV

echo "python extract_season.py $1 "  # $1 <-- season
python extract.py $1 > LOGS/extract_season_$1.log
