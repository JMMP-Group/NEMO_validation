#!/bin/bash
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --mem=20000
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --account=jmmp
module add jaspy
source activate $CONDA_ENV

python compare_surface_state.py
