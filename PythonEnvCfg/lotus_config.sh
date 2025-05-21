#!/bin/bash

## Set the machine to be used. Pick one.
export MACHINE="LOTUS"  # resource on JASMIN

## Years to loop over during monthly preprocessing: called in iter_en4_proc.sh, iter_sub_METEST.sh
export STARTYEAR=2004 #1980  # 2004
export ENDYEAR=2014 #1980    # 2014

## Process monthly data. Required in iter_sub_METEST.sh
export MOD="P2.0"  # Model reference name
export COMP_MOD="co7"

export MASS_DIR="xu-ct872"

# Use to activate conda environment in $MACHINE_pre_process_en4_monthly.sh
#export CONDA_ENV="/home/users/jelt/miniconda3/envs/coast_dev"
export CONDA_ENV="/home/users/ryapat30/.conda/envs/coast"

# location of COAsT repo, if using a particular branch
export COAST_REPO="/home/users/ryapat30/NOC/COAsT"
# COAsT configuration files
export DN_CFG="/home/users/ryapat30/NOC/COAsT/config/"
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
export DN_DOM="/gws/nopw/j04/jmmp/jmmp_collab/AMM15/DOMAIN_CFG/"

# directory for NEMO data files
#export DN_DAT="/gws/nopw/j04/jmmp/CO9_AMM15/outputs/p0/daily/"  # Dave 25h data
if [[ $MOD = "co7" ]]; then
 export DN_DAT="/gws/nopw/j04/jmmp/CO9_AMM15/outputs/co7/daily/"  # co7 24h ave
export GRID="/gws/nopw/j04/jmmp/public/AMM15/DOMAIN_CFG/CO7_EXACT_CFG_FILE.nc"  # contains the grid information for NEMO
else
 export DN_DAT="/gws/nopw/j04/jmmp/MASS/"$MASS_DIR"/daily/"  # MOD 24h ave
 if [[ $MOD = "P2.0" ]]; then
   export GRID="/gws/nopw/j04/jmmp/public/AMM15/DOMAIN_CFG/GEG_SF12.nc"
 else
   echo "WARNING: grid not set find matching grid and add to lotus_config.sh"
   # options?: "domain_cfg_MEs_01-003_opt_v1.nc" #  "GEG_SF12.nc"
   # options?: "GEG_SF12.nc" CO7_EXACT_CFG_FILE.nc
 fi
fi

# temporary fix for adding comparison model
if [[ $COMP_MOD = "co7" ]]; then
 export COMP_DAT="/gws/nopw/j04/jmmp/CO9_AMM15/outputs/co7/daily/"  # co7 24h ave
export COMP_GRID="/gws/nopw/j04/jmmp/public/AMM15/DOMAIN_CFG/CO7_EXACT_CFG_FILE.nc"  # contains the grid information for NEMO
else
 export COMP_DAT="/gws/nopw/j04/jmmp/MASS/"$MASS_DIR"/daily/"  # MOD 24h ave
 if [[ $COM_MOD = "P2.0" ]]; then
   export COMP_GRID="/gws/nopw/j04/jmmp/public/AMM15/DOMAIN_CFG/GEG_SF12.nc"
 else
   echo "WARNING: grid not set find matching grid and add to lotus_config.sh"
   # options?: "domain_cfg_MEs_01-003_opt_v1.nc" #  "GEG_SF12.nc"
   # options?: "GEG_SF12.nc" CO7_EXACT_CFG_FILE.nc
 fi
fi
# directory for analysis output
#export DN_OUT="/home/users/jelt/tmp/"$REGION"/" 
export DN_OUT="/gws/nopw/j04/jmmp/CO9_AMM15_validation/"$MOD"/"
export COMP_OUT="/gws/nopw/j04/jmmp/CO9_AMM15_validation/"$COMP_MOD"/"
