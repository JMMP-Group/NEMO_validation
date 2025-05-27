#!/bin/bash
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --account=jmmp
#SBATCH --mem=40000
#SBATCH -o LOGS/%A_%a.out
#SBATCH -e LOGS/%A_%a.err
#SBATCH --time=60:00
#SBATCH --ntasks=1
module add jaspy
source activate $CONDA_ENV

#time -v python  map_profiles.py $1 $2  $3 $4 > OUTPUTS/OUT_$1_$2_$3_$4
#time -v python  map_profiles.py $1 $2  $3 $4  > OUTPUTS/OUT_$1_$2_$3_$4
echo $1 $2 $3 $4 
python  map_profiles.py $1 $2 $3 $4  > LOGS/OUT_$1_$2_$3_$4_$5.log
