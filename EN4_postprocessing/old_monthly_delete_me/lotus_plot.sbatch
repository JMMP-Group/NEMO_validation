#!/bin/bash
#SBATCH --partition=short-serial
#SBATCH --mem=20000
#SBATCH --ntasks=1
#SBATCH -o LOGS/%A_%a_plot.out
#SBATCH -e LOGS/%A_%a_plot.err
#SBATCH --time=5
module add jaspy
source activate $CONDA_ENV

echo "python plot_month.py"
python plot_month.py > LOGS/plot_month.log