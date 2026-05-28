#!/bin/bash
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --account=jmmp
#SBATCH --mem=20000
#SBATCH -o LOGS/%A_%a.out
#SBATCH -e LOGS/%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
module add jaspy
source activate $CONDA_ENV

echo "python regional_masking.py " 
python regional_masking.py $1 > LOGS/regional_masking.log
