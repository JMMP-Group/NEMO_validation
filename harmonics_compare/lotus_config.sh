#export MACHINE="LOTUS"  # resource on JASMIN. Already set

## Years to loop over during preprocessing: called in ...
export STARTYEAR=2004 #1980  # 2004
export ENDYEAR=2014 #1980    # 2014

export CONDA_ENV="/home/users/jelt/miniconda3/envs/coast_dev"
# location of COAsT repo, if using a particular branch
export COAST_REPO="/home/users/jelt/GitHub/COAsT"

## RUN_NAME set in ${MACHINE,,}_pre_process_harm.sh from a passed variable
# E.g. ${MACHINE,,}_pre_process_harm.sh $config
#export RUN_NAME="GS1p1_tide"
#export RUN_NAME="GS1p2_full"
#export RUN_NAME="FES2014"


if [ $RUN_NAME = "GS1p1_tide" ]; then
 export FN_NEMO_DATA="/gws/nopw/j04/class_vol2/senemo/cwilso01/ZPS_REF_TIDE/OUTPUTS/SENEMO_1y_19810101_19811231_grid_T_2D.nc"

elif [ $RUN_NAME = "GS1p2_full" ]; then
 export FN_NEMO_DATA="/gws/nopw/j04/class_vol2/senemo/jdha/FINAL_TESTING/EXP_MESv2_NOTAPER_WAV_DJC_NTM_TDISSx2/output/SENEMO_1y_19810101_19811231_grid_T_2D.nc"

elif [ $RUN_NAME = "FES2014" ]; then
  export FN_FES_AMP="/gws/nopw/j04/class_vol2/senemo/shared/FES2014_M2_amp_on_SENEMO_grid.nc"
  export FN_FES_PHA="/gws/nopw/j04/class_vol2/senemo/shared/FES2014_M2_pha_on_SENEMO_grid.nc"
#else #

fi



export FN_HARM_OBS="/home/users/jelt/data/obs/for_DA_dense/obs_M2.nc"

export FN_NEMO_DOMAIN="/gws/nopw/j04/class_vol2/senemo/jdha/FINAL_TESTING/EXP_MESv2_NOTAPER_WAV_DJC_NTM_TDISSx2/config/domain_cfg.nc"
export FN_NEMO_CFG="/home/users/jelt/GitHub/COAsT/config/example_nemo_grid_t.json"

export DN_OUT="/gws/nopw/j04/class_vol2/senemo/jelt/"

export FN_ANALYSIS_OUT=$DN_OUT"PROCESSED/"$RUN_NAME".nc"



