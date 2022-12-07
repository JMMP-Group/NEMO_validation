#! /bin/ksh
#SBATCH --mem=20000
#SBATCH --ntasks=1
#SBATCH --output=mean_output
#SBATCH --error=mean_error
#SBATCH --time=5

time -v ipython mean_monthly.py $1 $2 #(Model month)


