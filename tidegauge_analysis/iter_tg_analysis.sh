#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
source config.sh

rm LOGS/OUT* LOGS/*.err LOGS/*.out

for (( port=0; port<$(expr $N_PORTS + 1); port++ ))
do
  echo "sbatch -J id_dim"$port" "${MACHINE,,}"_rs_hourly_ssh.pbs" $port
  sbatch -J id_dim${port} ${MACHINE,,}_rs_hourly_ssh.pbs $port
done
