#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
for (( start=1980; start<1981; start++ ))
do
  for (( month=1; month<13; month++ ))
  do 
     end=$(expr $start + 1 )
#     sbatch spice_ana_SEAS_P0.0.sh $SEAS $start $end
     #echo "sbatch spice_ana_MOD.sh $start $month"
     #sbatch spice_pre_process_en4_monthly.sh $start $month
     echo "sbatch lotus_pre_process_en4_monthly.sh $start $month"
     sbatch lotus_pre_process_en4_monthly.sh $start $month
  done
done
