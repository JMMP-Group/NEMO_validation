#!/bin/bash 
#SBATCH --partition=short-serial 
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --time=30:00
#SBATCH --array=1-132%10
module add jaspy
source activate /home/users/jelt/miniconda3/envs/coast_dev
#source /home/users/jelt/miniconda3/envs/coast_dev/lib/python3.9/venv/scripts/common/activate
#source ~/envs/coast/bin/activate
# Pass in integer year and month: python pre_process_en4_monthly.py <year> <month>
python pre_process_en4_monthly.py $1 $2
