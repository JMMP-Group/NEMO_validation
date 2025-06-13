#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
source config.sh  # sets $MACHINE. Is called again before other variables are used

rm LOGS/OUT* LOGS/*.err LOGS/*.out

echo "sbatch "${MACHINE,,}"_pre_process_obs.sh"
sbatch -J "pre_proc" ${MACHINE,,}_pre_process_obs.sh --exclude=host1185

