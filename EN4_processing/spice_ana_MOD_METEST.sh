#!/bin/bash -l
#SBATCH --mem=40000
#SBATCH --ntasks=1
#SBATCH --output=ana_output
#SBATCH --error=ana_error
#SBATCH --time=180
source config.sh
conda activate $CONDA_ENV 

#time -v python  GEN_MOD_Dave_example_profile_validation_SEAS.py $1 $2  $3 $4 > OUTPUTS/OUT_$1_$2_$3_$4
#time -v python  GEN_MOD_Dave_example_profile_validation_METEST.py $1 $2  $3 $4  > OUTPUTS/OUT_$1_$2_$3_$4
echo $1 $2 $3 $4 $4
time  python -v GEN_MOD_Dave_example_profile_validation.py $1 $2 $3 $4 $5  > LOGS/OUT_$1_$2_$3_$4_$5.log
