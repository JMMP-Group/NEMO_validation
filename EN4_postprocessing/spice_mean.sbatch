#! /bin/ksh
#SBATCH --mem=20000
#SBATCH --ntasks=1
#SBATCH -o LOGS/%A_%a_mean.out
#SBATCH -e LOGS/%A_%a_mean.err
#SBATCH --time=5

echo "time -v python mean_monthly.py $1 $2"  # $1 <-- MOD, $2 <-- month
# NB was ipython
time -v python mean_monthly.py $1 $2 > LOGS/mean_monthly_$1_$2.log
