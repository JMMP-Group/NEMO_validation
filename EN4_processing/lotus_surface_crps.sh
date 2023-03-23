#!/bin/bash
#SBATCH --partition=short-serial
#SBATCH --mem=40000
#SBATCH -o LOGS/%A_%a.out
#SBATCH -e LOGS/%A_%a.err
#SBATCH --time=60:00
#SBATCH --ntasks=1
module add jaspy
source activate $CONDA_ENV


echo $1 $2 $3 $4 $5
python  surface_crps.py $1 $2 $3 $4 $5  > LOGS/OUT_$1_$2_$3_$4_$5.log
