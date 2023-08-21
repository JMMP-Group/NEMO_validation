#!/bin/bash
#SBATCH --partition=short-serial
#SBATCH --mem=20000
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
module add jaspy
source activate $CONDA_ENV

python compare_surface_state.py
