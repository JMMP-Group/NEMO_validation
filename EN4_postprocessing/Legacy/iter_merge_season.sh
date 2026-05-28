#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
source config.sh

mkdir -p LOGS
rm LOGS/*err LOGS/*out LOGS/*log

for season in DJF MAM JJA SON;
do
     echo "sbatch -J "${MOD}${season} ${MACHINE,,}"_merge_season.sbatch $MOD $season"
     sbatch -J ${MOD}${season} ${MACHINE,,}_merge_season.sbatch $MOD $season
done
