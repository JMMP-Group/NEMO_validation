#! /bin/ksh
#SBATCH --mem=20000
#SBATCH --ntasks=1
#SBATCH --output=ana_output
#SBATCH --error=ana_error
#SBATCH --time=5

time -v python merge_monthly.py $1 $2 #(Model month)


