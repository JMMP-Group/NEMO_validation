#export MACHINE="LOTUS"  # resource on JASMIN. Already set

## Years to loop over during monthly preprocessing: called in iter_tg_analysis.sh
export STARTYEAR=2004 #1980  # 2004
export ENDYEAR=2014 #1980    # 2014

export CONDA_ENV="/home/users/jelt/miniconda3/envs/coast_dev"
# location of COAsT repo, if using a particular branch
export COAST_REPO="/home/users/jelt/GitHub/COAsT"


export RUN_NAME="P0.0"
export RUN_NAME="P1.5b"
export RUN_NAME="P1.5"
export RUN_NAME="P1.5c"
export RUN_NAME="co7"
export RUN_NAME="P2.0"

## 2014 P0.0
#export MASS="xu-cb676"
#export FN_NEMO_DATA="/gws/nopw/j04/jmmp/MASS/"$MASS"/shelftmb/YYYYMM*_shelftmb_grid_T.nc"

# 10 years P0.0
if [ $RUN_NAME = "P0.0" ]; then
 export FN_NEMO_DATA="/gws/nopw/j04/jmmp_collab/AMM15/OUTPUTS/"$RUN_NAME"/NC3_SSH/YYYYMM*_shelftmb_grid_T.nc"

elif [ $RUN_NAME = "co7" ]; then
 export FN_NEMO_DATA="/gws/nopw/j04/jmmp/MASS/rosie_mi-ao113_2004/field.nc.file/YYYYMM*_shelftmb_grid_T.nc"

else # elif [ $RUN_NAME = "P1.5" ]; then
 export FN_NEMO_DATA="/gws/nopw/j04/jmmp_collab/AMM15/OUTPUTS/"$RUN_NAME"/SSH/YYYYMM*_shelftmb_grid_T.nc.ppc3"
fi



export FN_NEMO_DOMAIN="/gws/nopw/j04/jmmp/jmmp_collab/AMM15/DOMAIN_CFG/CO7_EXACT_CFG_FILE.nc"
export FN_NEMO_CFG="/home/users/jelt/GitHub/COAsT/config/example_nemo_grid_t.json"

export FN_OUT="/gws/nopw/j04/jmmp/CO9_AMM15_validation/"$RUN_NAME"/tg_analysis/ssh_hourly_"$RUN_NAME".nc"
export FN_OBS="/gws/nopw/j04/jmmp/CO9_AMM15/obs/tg_amm15.nc"
export N_PORTS=61 # number of ports in $FN_OBS

export DN_OUT="/gws/nopw/j04/jmmp/CO9_AMM15_validation/"$RUN_NAME"/tg_analysis/FIGS"

## Use validate_ssh_tg_hourly.py
export FN_EXTR_OUT="/gws/nopw/j04/jmmp/CO9_AMM15_validation/"$RUN_NAME"/tg_analysis/ssh_hourly_extract_"$RUN_NAME".nc"
export FN_ANALYSE_OUT="/gws/nopw/j04/jmmp/CO9_AMM15_validation/"$RUN_NAME"/tg_analysis/ssh_hourly_analyse_"$RUN_NAME".nc"
export FN_SSH_HOURLY_STATS="/gws/nopw/j04/jmmp/CO9_AMM15_validation/"$RUN_NAME"/tg_analysis/ssh_hourly_analyse_"$RUN_NAME".nc"
