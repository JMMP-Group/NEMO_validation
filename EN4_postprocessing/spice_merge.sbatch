#!/bin/ksh
#SBATCH --mem=20000
#SBATCH --ntasks=1
#SBATCH --output=LOGS/ana_output
#SBATCH --error=LOGS/ana_error
#SBATCH --time=5

echo "time -v python merge_monthly.py $1 $2"  # $1 <-- MOD, $2 <-- month
time -v python merge_monthly.py $1 $2 > LOGS/merge_monthly_$1_$2.log