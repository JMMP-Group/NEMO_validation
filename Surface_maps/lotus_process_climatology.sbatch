#!/bin/bash
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --mem=500G
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --account=jmmp
source activate nemo_validation

python -u compare_surface_state.py
