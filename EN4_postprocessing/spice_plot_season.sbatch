#!/bin/bash -l
#SBATCH --mem=20000
#SBATCH --ntasks=1
#SBATCH --output LOGS/%A_%a_plot.out
#SBATCH --error LOGS/%A_%a_plot.err
#SBATCH --time=5
source config.sh
conda activate $CONDA_ENV 

echo "time python plot_season.py"
time python plot_season.py > LOGS/plot_season.log