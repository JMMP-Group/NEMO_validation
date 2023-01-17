#!/bin/bash

## Set the machine to be used. Pick one.
export MACHINE="LOTUS"  # resource on JASMIN

## Years to loop over during monthly preprocessing: called in iter_en4_proc.sh, iter_sub_METEST.sh
export STARTYEAR=2004 #1980  # 2004
export ENDYEAR=2014 #1980    # 2014

## Process monthly data. Required in iter_sub_METEST.sh
export MOD="P0.0"  # Model reference name
export GRID="CO7_EXACT_CFG_FILE.nc"  # contains the grid information for NEMO
# options?: "domain_cfg_MEs_01-003_opt_v1.nc" #  "GEG_SF12.nc"
# options?: "GEG_SF12.nc" CO7_EXACT_CFG_FILE.nc
if [[ $MOD = "P0.0" ]]
then
 export MASS_DIR="xu-cb676"
else
 export MASS_DIR=""
 echo "Not ready for that experiment"
fi

# Use to activate conda environment in $MACHINE_pre_process_en4_monthly.sh
export CONDA_ENV="/home/users/jelt/miniconda3/envs/coast_dev"

# location of COAsT repo, if using a particular branch
export COAST_REPO="/home/users/jelt/GitHub/COAsT"
# COAsT configuration files
export DN_CFG="/home/users/jelt/GitHub/COAsT/config/"
export FN_CFG_PROF=$DN_CFG"example_en4_profiles.json"
export FN_CFG_NEMO=$DN_CFG"example_nemo_grid_t.json"

# location of raw EN4 data
#export DIN_EN4="/home/users/jelt/EN4/"
export DIN_EN4="/gws/nopw/j04/class_vol2/senemo/shared/EN4/downloads/EN.4.2.2.profiles/"

# prefix for preprocessed EN4 data (chunked into files by region and time)
export REGION="AMM15"
# location of preprocessed EN4 data
export DOUT_EN4="/gws/nopw/j04/class_vol2/senemo/shared/EN4/processed/"$REGION"/"

# directory for NEMO domain_cfg.nc file
export DN_DOM="/gws/nopw/j04/jmmp_collab/CO9_AMM15/inputs/domains/"
# directory for NEMO data files
#export DN_DAT="/gws/nopw/j04/jmmp/CO9_AMM15/outputs/p0/daily/"  # Dave 25h data
export DN_DAT="/gws/nopw/j04/jmmp/MASS/"$MASS_DIR"/daily/"  # P0.0 24h ave

# directory for analysis output
#export DN_OUT="/home/users/jelt/tmp/"$REGION"/" 
export DN_OUT="/gws/nopw/j04/jmmp/CO9_AMM15_validation/"$MOD"/profiles/"
