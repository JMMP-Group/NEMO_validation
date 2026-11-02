#export MACHINE="LOTUS"  # resource on JASMIN. Already set

## Years to loop over during monthly preprocessing: called in iter_tg_analysis.sh
export STARTYEAR=2013 #1980  # 2004
export ENDYEAR=2013 #1980    # 2014

export CONDA_ENV="/home/users/o.lambertbrown/.conda/envs/nemo_valid"  ## OR WHAT IS THIS REALLY? WHAT IS THE PATH?
# location of COAsT repo, if using a particular branch
export COAST_REPO="/home/users/o.lambertbrown/COAsT"



export RUN_NAME="P1.5c"

## 2014 P0.0
#export MASS="xu-cb676"
#export FN_NEMO_DATA="/gws/nopw/j04/jmmp/MASS/"$MASS"/shelftmb/YYYYMM*_shelftmb_grid_T.nc"

export FN_NEMO_DATA='/data/users/o.lambertbrown/Datasets/models/CO9/u-cu674/YYYYMM*_shelftmb_grid_T.nc.ppc3'


export DN_CFG="/home/users/o.lambertbrown/COAsT/config/"
export FN_NEMO_CFG=$DN_CFG"example_nemo_grid_t.json"

export FN_NEMO_DOMAIN="/data/users/o.lambertbrown/Datasets/models/CO9/u-cu674/domain_cfg_sf12.nc"

export FN_OUT="/data/users/o.lambertbrown/Datasets/models/CO9/u-cu674/analysis/tides/ssh_hourly_"$RUN_NAME".nc"
export FN_OBS="/data/scratch/o.lambertbrown/GESLA/AMM15_15m_tidegauge_merged.nc"
export N_PORTS=153 # number of ports in $FN_OBS

export DN_OUT="/home/users/o.lambertbrown/.conda/envs/coast-test/NEMO_validation/Tidegauge_processing/FIGS"

## Use validate_ssh_tg_hourly.py
export FN_EXTR_OUT="/data/users/o.lambertbrown/Datasets/models/CO9/u-cu674/analysis/tides/ssh_hourly_extract_"$RUN_NAME".nc"
export FN_ANALYSE_OUT="/data/users/o.lambertbrown/Datasets/models/CO9/u-cu674/analysis/tides/ssh_hourly_analyse_"$RUN_NAME".nc"
export FN_SSH_HOURLY_STATS="/data/users/o.lambertbrown/Datasets/models/CO9/u-cu674/analysis/tides/ssh_hourly_corr_"$RUN_NAME".nc"
