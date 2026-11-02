#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
source config.sh

rm LOGS/OUT* LOGS/*.err LOGS/*.out

for (( port=0; port<$(expr $N_PORTS + 1); port++ ))
do
    echo "sbatch -J port"${port}" "${MACHINE,,}"_process_en4_monthly.sh "$port
    sbatch -J port${port} ${MACHINE,,}_process_tg_monthly.sh $port
done
