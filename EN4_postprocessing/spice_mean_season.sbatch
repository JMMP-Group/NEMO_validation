#! /bin/ksh
#SBATCH --mem=20000
#SBATCH --ntasks=1
#SBATCH -o LOGS/%A_%a_mean.out
#SBATCH -e LOGS/%A_%a_mean.err
#SBATCH --time=30:00

echo "time -v python mean_season.py $1 $2"  # $1 <-- MOD, $2 <-- season
# NB was ipython
time -v python mean_season.py $1 $2 > LOGS/mean_season_$1_$2.log
