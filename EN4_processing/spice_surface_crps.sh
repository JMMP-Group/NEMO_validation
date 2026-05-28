#!/bin/bash -l
#SBATCH --mem=40000
#SBATCH --output=LOGS/ana_crps$1$2.out
#SBATCH --error=LOGS/ana_crps$1$2.err
#SBATCH --time=180
#SBATCH --ntasks=1
source config.sh
conda activate $CONDA_ENV

start=$(date +%s)

echo $1 $2 
time python -v extract_model_surface.py $1 $2 > LOGS/OUT_$1_$2.log

mid=$(date +%s)
echo "Extracted. Time taken: $(($mid - $start)) s" > LOGS/$1_$2_timing.log

time python -v surface_crps.py $1 $2 >> LOGS/OUT_$1_$2.log

end=$(date +%s)
echo "crps processing time: $(($end - $start)) s" >> LOGS/$1_$2_timing.log
