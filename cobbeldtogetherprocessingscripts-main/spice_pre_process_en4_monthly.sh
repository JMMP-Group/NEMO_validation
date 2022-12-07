#! /bin/ksh
#SBATCH --mem=200000
#SBATCH --ntasks=1
#SBATCH --output=pre_en4_output
#SBATCH --error=pre_en4_error
#SBATCH --time=180
source config.sh
source activate $CONDA_ENV

#python pre_process_en4_1990_2020.py
# Pass in integer year and month: python pre_process_en4_monthly.py <year> <month>
python pre_process_en4_monthly.py $1 $2
