#!/bin/ksh
#SBATCH --mem=20000
#SBATCH --ntasks=1
#SBATCH -o LOGS/%A_%a_plot.out
#SBATCH -e LOGS/%A_%a_plot.err
#SBATCH --time=5

echo "time -v python plot_monthy.py"
time -v python plot_month.py > LOGS/plot_month.log