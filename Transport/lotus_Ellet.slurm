#!/bin/bash

##BATCH --account=jmmp
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --account=jmmp
#SBATCH --mem=50GB
#SBATCH --ntasks=1
#SBATCH --time=08:00:00

#module load jaspy
conda activate nemo_validation

python -u plot_transport_section.py 

