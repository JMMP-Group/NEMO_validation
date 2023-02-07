#export MACHINE="LOTUS"  # resource on JASMIN. Already set

export CONDA_ENV="/home/users/jelt/miniconda3/envs/coast_dev"
# location of COAsT repo, if using a particular branch
export COAST_REPO="/home/users/jelt/GitHub/COAsT"

export RUN_NAME="P0.0"
export RUN_NAME="P0.0_test"
export MASS="xu-cb676"

# 2014 P0.0
export FN_NEMO_DATA="/gws/nopw/j04/jmmp/MASS/"$MASS"/shelftmb/2014*_shelftmb_grid_T.nc"
# 10 years P0.0
export FN_NEMO_DATA="/gws/nopw/j04/jmmp_collab/AMM15/OUTPUTS/P0.0/NC3_SSH/20*_shelftmb_grid_T.nc"

export FN_NEMO_DOMAIN="/gws/nopw/j04/jmmp_collab/CO9_AMM15/inputs/domains/CO7_EXACT_CFG_FILE.nc"
export FN_NEMO_CFG="/home/users/jelt/GitHub/COAsT/config/example_nemo_grid_t.json"

export FN_OUT="/gws/nopw/j04/jmmp/CO9_AMM15_validation/"$RUN_NAME"/tg_analysis/ssh_hourly_"$RUN_NAME".nc"
export FN_OBS="/gws/nopw/j04/jmmp/CO9_AMM15/obs/tg_amm15.nc"

export DN_OUT="/home/users/jelt/GitHub/NEMO_validation/tidegauge_analysis/FIGS"

## Use validate_ssh_tg_hourly.py
export FN_EXTR_OUT="/gws/nopw/j04/jmmp/CO9_AMM15_validation/"$RUN_NAME"/tg_analysis/ssh_hourly_extract_"$RUN_NAME".nc"
export FN_ANALYSE_OUT="/gws/nopw/j04/jmmp/CO9_AMM15_validation/"$RUN_NAME"/tg_analysis/ssh_hourly_analyse_"$RUN_NAME".nc"
export FN_SSH_HOURLY_STATS="/gws/nopw/j04/jmmp/CO9_AMM15_validation/"$RUN_NAME"/tg_analysis/ssh_hourly_analyse_"$RUN_NAME".nc"
