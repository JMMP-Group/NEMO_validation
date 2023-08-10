#!/bin/bash
#SBATCH --partition=short-serial
#SBATCH --mem=40000
#SBATCH -o LOGS/%A_%a.out
#SBATCH -e LOGS/%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
module add jaspy
source activate $CONDA_ENV

start=$(date +%s)

echo $1 $2 
python extract_model_surface.py $1 $2 > LOGS/OUT_$1_$2.log

mid=$(date +%s)
echo "Extracted. Time taken: $(($mid - $start)) s" > LOGS/$1_$2_timing.log

python surface_crps.py $1 $2 >> LOGS/OUT_$1_$2.log

end=$(date +%s)
echo "crps processing time: $(($end - $start)) s" >> LOGS/$1_$2_timing.log
