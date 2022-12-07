#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
source config.sh

#for (( start=1980; start<1981; start++ ))
for (( start=$STARTYEAR; start<$(expr $ENDYEAR + 1); start++ ))
do
  for (( month=1; month<13; month++ ))
  do
    if [ $MACHINE = "SPICE" ]; then
      echo "sbatch "+${MACHINE,,}+"_pre_process_en4_monthly.sh $start $month"
      sbatch spice_pre_process_en4_monthly.sh $start $month
    elif [ $MACHINE = "LOTUS" ]; then
      echo "sbatch "+${MACHINE,,}+"_pre_process_en4_monthly.sh $start $month"
      sbatch ${MACHINE,,}_pre_process_en4_monthly.sh $start $month
    else
      echo "not expecting $MACHINE"
    fi
  done
done
