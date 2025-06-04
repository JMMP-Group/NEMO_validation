#!/bin/bash
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --account=jmmp
#SBATCH --mem=20000
#SBATCH --ntasks=1
#SBATCH -o LOGS/%A_%a_plot.out
#SBATCH -e LOGS/%A_%a_plot.err
#SBATCH --time=20:00
module add jaspy
conda activate $CONDA_ENV

echo "python plot_season.py"
python plot_season.py > LOGS/plot_season.log
