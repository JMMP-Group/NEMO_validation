#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
source config.sh

rm LOGS/OUT* LOGS/*.err LOGS/*.out

#MOD=$1
#GRID="domain_cfg_MEs_01-003_opt_v1.nc" #  "GEG_SF12.nc"
#GRID=$2 #   "GEG_SF12.nc" CO7_EXACT_CFG_FILE.nc
#for (( start=2004; start<2015; start++ ))
for (( start=$STARTYEAR; start<$(expr $ENDYEAR + 1); start++ ))

  do
  for (( month=1; month<13; month++ ))
  #for (( month=1; month<2; month++ ))
    do 
     end=$(expr $start + 1 )
#     sbatch spice_ana_SEAS_P0.0.sh $SEAS $start $end
     #sbatch spice_ana_MOD_METEST.sh $MOD $start $end $GRID
     #sbatch  spice_ana_MOD_METEST.sh $MOD $start $month $end $GRID
     echo "sbatch "${MACHINE,,}"_ana_MOD_METEST.sh $MOD $start $month $end $GRID"
     sbatch ${MACHINE,,}_ana_MOD_METEST.sh $MOD $start $month $end $GRID
#     sh spice_ana_MOD_METEST.sh $MOD $start $end $GRID
 done
done
