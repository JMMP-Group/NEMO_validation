#!/bin/bash 
#SBATCH --partition=short-serial 
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --time=30:00
#SBATCH --array=1-132%10

## THIS FILES NEEDS EDITTING TO RUN ON SPICE AT THE MO

module add jaspy
source activate $CONDA_ENV

# Pass in integer year and month: python pre_process_en4_monthly.py <year> <month>
python pre_process_en4_monthly.py $1 $2
