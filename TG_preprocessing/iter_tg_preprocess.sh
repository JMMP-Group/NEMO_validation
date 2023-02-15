#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
source config.sh

rm LOGS/OUT* LOGS/*.err LOGS/*.out

#for (( start=2004; start<2005; start++ ))
for (( start=$STARTYEAR; start<$(expr $ENDYEAR + 1); start++ ))
do
  for (( month=1; month<13; month++ ))
  do
    echo "sbatch "${MACHINE,,}"_pre_process_en4_monthly.sh $start $month"
    sbatch -J ${start}${month} ${MACHINE,,}_pre_process_tg_monthly.sh $start $month
  done
done
