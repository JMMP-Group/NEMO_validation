#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
cd ../PythonEnvCfg/
source config.sh
cd ../EN4_postprocessing
conda activate $CONDA_ENV

mkdir -p LOGS
rm LOGS/*err LOGS/*out LOGS/*log

for season in DJF MAM JJA SON;
do
   echo "sbatch -J "${MOD}${season} ${MACHINE,,}"_extract_season.sbatch $season"
   sbatch -J ${MOD}${season} ${MACHINE,,}_extract_season.sbatch $season
done
