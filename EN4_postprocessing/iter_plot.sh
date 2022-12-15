#!/bin/bash
#SBATCH --partition=short-serial
#SBATCH --mem=20000
#SBATCH --ntasks=1
#SBATCH -o LOGS/%A_%a_plot.out
#SBATCH -e LOGS/%A_%a_plot.err
#SBATCH --time=5

echo "Bash version ${BASH_VERSION}..."
source config.sh

mkdir -p LOGS
rm LOGS/*err LOGS/*out LOGS/*log

mkdir -p FIGS

module add jaspy
source activate $CONDA_ENV

#for (( month=1; month<13; month++ ))
#do
#     echo "sbatch "${MACHINE,,}"_mean.sbatch $MOD $month"
#     sbatch ${MACHINE,,}_mean.sbatch $MOD $month
#done
echo "plot_month.py"
python plot_month.py > LOGS/plot_month.log