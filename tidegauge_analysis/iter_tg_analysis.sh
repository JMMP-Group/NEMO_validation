#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
source config.sh

rm LOGS/OUT* LOGS/*.err LOGS/*.out

echo "sbatch "${MACHINE,,}"_rs_hourly_ssh.pbs"
sbatch ${MACHINE,,}_rs_hourly_ssh.pbs