#!/bin/bash

##BATCH --account=jmmp
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --account=jmmp
#SBATCH --mem=200G
#SBATCH --ntasks=1
#SBATCH --time=24:00:00

#module load jaspy
conda activate nemo_validation

python -u plot_transport_section.py 

