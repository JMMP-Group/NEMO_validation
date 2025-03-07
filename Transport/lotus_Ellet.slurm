#!/bin/bash

##BATCH --account=jmmp
#SBATCH --partition=short-serial
#SBATCH --mem=50GB
#SBATCH --ntasks=1
#SBATCH --time=03:00:00
##BATCH --qos=short
module add jaspy
source activate nemo_validation

python -u plot_transport_section.py 

