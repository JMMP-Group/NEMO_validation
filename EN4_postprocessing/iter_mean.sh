#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
source config.sh

mkdir -p LOGS
rm LOGS/*err LOGS/*out LOGS/*log

#for (( month=1; month<13; month++ ))
for (( month=1; month<2; month++ ))
do 
     echo "sbatch "${MACHINE,,}"_mean.sbatch $MOD $month"
     sbatch ${MACHINE,,}_mean.sbatch $MOD $month
done