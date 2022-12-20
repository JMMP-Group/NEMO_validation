#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
source config.sh

mkdir -p LOGS
rm LOGS/*err LOGS/*out LOGS/*log

mkdir -p FIGS

echo "sbatch "${MACHINE,,}"_plot_season.sbatch"
sbatch ${MACHINE,,}_plot_season.sbatch


