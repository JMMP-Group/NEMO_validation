#!/bin/bash -l
echo "Bash version ${BASH_VERSION}..."
source config.sh
conda activate $CONDA_ENV

mkdir -p LOGS
rm LOGS/*err LOGS/*out LOGS/*log

for season in DJF MAM JJA SON;
do
     echo "sbatch -J "${MOD}${season} ${MACHINE,,}"_extract_season.sbatch $MOD $season"
     sbatch -J ${MOD}${season} ${MACHINE,,}_extract_season.sbatch $MOD $season
done
