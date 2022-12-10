#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
source config.sh

mkdir -p LOGS

#MOD=$1
for (( month=1; month<13; month++ ))
do 
     echo "sbatch "${MACHINE,,}"_merge.sbatch $MOD $month"
     sbatch ${MACHINE,,}_merge.sbatch $MOD $month
done
