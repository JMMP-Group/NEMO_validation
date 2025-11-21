#!/bin/bash
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --mem=100GB
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --account=jmmp
source activate nemo_validation

python compare_surface_state.py
