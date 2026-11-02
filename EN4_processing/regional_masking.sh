#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
cd ../PythonEnvCfg/
source config.sh
cd ../EN4_processing
conda activate $CONDA_ENV

mkdir -p LOGS
rm LOGS/*err LOGS/*out LOGS/*log

# create regional masking
sbatch -J ${MOD}_regional_mask ${MACHINE,,}_regional_masking.slurm 
