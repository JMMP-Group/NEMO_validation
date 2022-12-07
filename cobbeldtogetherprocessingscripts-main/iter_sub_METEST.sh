
#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
MOD=$1
#GRID="domain_cfg_MEs_01-003_opt_v1.nc" #  "GEG_SF12.nc"
GRID=$2 #   "GEG_SF12.nc" CO7_EXACT_CFG_FILE.nc
 for (( start=2004; start<2015; start++ ))
  do
  for (( month=1; month<13; month++ ))
  #for (( month=12; month<13; month++ ))
    do 
     end=$(expr $start + 1 )
#     sbatch spice_ana_SEAS_P0.0.sh $SEAS $start $end
     #sbatch spice_ana_MOD_METEST.sh $MOD $start $end $GRID
     sbatch  spice_ana_MOD_METEST.sh $MOD $start $month $end $GRID
#     sh spice_ana_MOD_METEST.sh $MOD $start $end $GRID
 done
done
