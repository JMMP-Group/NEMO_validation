#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
source config.sh  # sets $MACHINE. Is called again before other variables are used

rm LOGS/OUT* LOGS/*.err LOGS/*.out

#configArray=("GS1p1_tide" "GS1p2_full" "FES2014")
configArray=("GS1p1_tide_final" "GS1p7_triads" "FES2014")
#configArray=("FES2014")
for config in ${configArray[*]}
do
  echo "sbatch "${MACHINE,,}"_pre_process_harm.sh $config"
  sbatch -J ${config} ${MACHINE,,}_pre_process_harm.sh $config --exclude=host1185
done
