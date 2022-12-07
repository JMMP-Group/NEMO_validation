#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
source config.sh

#for (( start=1980; start<1981; start++ ))
for (( start=$STARTYEAR; start<$(expr $ENDYEAR + 1); start++ ))
do
  for (( month=1; month<13; month++ ))
  do
     end=$(expr $start + 1 )
#     sbatch spice_ana_SEAS_P0.0.sh $SEAS $start $end
     if [ $MACHINE = "SPICE" ]; then
      echo "sbatch spice_pre_process_en4_monthly.sh $start $month"
      sbatch spice_pre_process_en4_monthly.sh $start $month  # THIS FILE IS MISSING
     elif [ $MACHINE = "LOTUS" ]; then
      echo "sbatch lotus_pre_process_en4_monthly.sh $start $month"
      sbatch lotus_pre_process_en4_monthly.sh $start $month
    fi
  done
done
