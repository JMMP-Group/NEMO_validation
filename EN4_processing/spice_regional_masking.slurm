#!/bin/bash
#SBATCH --mem=40000
#SBATCH -o LOGS/%A_%a.out
#SBATCH -e LOGS/%A_%a.err
#SBATCH --time=180
#SBATCH --ntasks=1
source config.sh
source activate $CONDA_ENV

echo "python regional_masking.py " 
python regional_masking.py $1 > LOGS/regional_masking.log
