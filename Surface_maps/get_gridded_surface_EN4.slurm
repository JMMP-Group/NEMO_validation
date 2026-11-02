#!/bin/bash
#SBATCH --partition=short-serial
#SBATCH --mem=64GB
#SBATCH -o LOGS/%A_%a.out
#SBATCH -e LOGS/%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
module add jaspy
source activate $CONDA_ENV
export PYTHONPATH=$PYTHONPATH:/home/users/ryapat30/NOC

start=$(date +%s)

python get_gridded_surface_EN4.py 2004 1 2015

end=$(date +%s)
echo "get_gridded_surface_EN4 processing time: $(($end - $start)) s" >> LOGS/$1_$2_timing.log
