#!/bin/bash
#SBATCH --partition=short-serial
######SBATCH --partition=high-mem
#SBATCH --mem=40000
#SBATCH -o LOGS/%A_%a.out
#SBATCH -e LOGS/%A_%a.err
#SBATCH --time=30:00
#SBATCH --ntasks=1
module add jaspy
source activate $CONDA_ENV

# SPICE
# X##! /bin/ksh
#X##SBATCH --mem=40000
#X##SBATCH --ntasks=1
#X##SBATCH --output=ana_output
#X##SBATCH --error=ana_error
#X##SBATCH --time=180

#time -v python  GEN_MOD_Dave_example_profile_validation_SEAS.py $1 $2  $3 $4 > OUTPUTS/OUT_$1_$2_$3_$4
#time -v python  GEN_MOD_Dave_example_profile_validation_METEST.py $1 $2  $3 $4  > OUTPUTS/OUT_$1_$2_$3_$4
echo $1 $2 $3 $4 $5
python  GEN_MOD_Dave_example_profile_validation.py $1 $2 $3 $4 $5  > LOGS/OUT_$1_$2_$3_$4_$5.log
