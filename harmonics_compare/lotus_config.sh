#export MACHINE="LOTUS"  # resource on JASMIN. Already set

export CONDA_ENV="/home/users/jelt/miniconda3/envs/coast_dev"
# location of COAsT repo, if using a particular branch
export COAST_REPO="/home/users/jelt/GitHub/COAsT"

## RUN_NAME set in ${MACHINE,,}_pre_process_harm.sh from a passed variable
# E.g. ${MACHINE,,}_pre_process_harm.sh $config
#export RUN_NAME="GS1p1_tide"
#export RUN_NAME="GS1p2_full"
#export RUN_NAME="FES2014"


if [ "$RUN_NAME" = "GS1p1_tide" ]; then
 export FN_NEMO_DATA="/gws/nopw/j04/class_vol2/senemo/RUNS2024r01/GS1p1_tide/output/SENEMO_1y_19810101_19811231_grid_T_2D.nc"
 export FN_NEMO_DOMAIN="/gws/nopw/j04/class_vol2/senemo/RUNS2024r01/GS1p1_tide/config/domain_cfg.nc"

elif [ "$RUN_NAME" = "GS1p2_full" ]; then
 export FN_NEMO_DATA="/gws/nopw/j04/class_vol2/senemo/RUNS2024r01/GS1p2_full/output/SENEMO_1y_19810101_19811231_grid_T_2D.nc"
 export FN_NEMO_DOMAIN="/gws/nopw/j04/class_vol2/senemo/RUNS2024r01/GS1p2_full/config/domain_cfg.nc"

elif [ "$RUN_NAME" = "GS1p6_full" ]; then
 export FN_NEMO_DATA="/gws/nopw/j04/class_vol2/senemo/jholt/EXP_GS1p6_full_IWD_soenhance_SSR_ice/SENEMO_1y_19810101_19811231_grid_T_2D.nc"
 export FN_NEMO_DOMAIN="/gws/nopw/j04/class_vol2/senemo/RUNS2024r01/GS1p2_full/config/domain_cfg.nc"

elif [ "$RUN_NAME" = "GS1p6_full_so" ]; then
 export FN_NEMO_DATA="/gws/nopw/j04/class_vol2/senemo/jholt/EXP_GS1p6_full_IWD_soenhance_SSR_ice_so/SENEMO_1y_19810101_19811231_grid_T_2D.nc"
 export FN_NEMO_DOMAIN="/gws/nopw/j04/class_vol2/senemo/RUNS2024r01/GS1p2_full/config/domain_cfg.nc"

elif [ "$RUN_NAME" = "GS1p7" ]; then
 export FN_NEMO_DATA="/gws/nopw/j04/class_vol2/senemo/jholt/EXP_GS1p6_full_IWD_soenhance_qdr2/SENEMO_1y_19810101_19811231_grid_T_2D.nc"
 export FN_NEMO_DOMAIN="/gws/nopw/j04/class_vol2/senemo/RUNS2024r01/GS1p2_full/config/domain_cfg.nc"

elif [ "$RUN_NAME" = "FES2014" ]; then
  export DN_FES="/gws/nopw/j04/class_vol2/senemo/shared/FES2014/"
  #export FN_FES_AMP="/gws/nopw/j04/class_vol2/senemo/shared/FES2014_M2_amp_on_SENEMO_grid.nc"
  #export FN_FES_PHA="/gws/nopw/j04/class_vol2/senemo/shared/FES2014_M2_pha_on_SENEMO_grid.nc"
#else #

elif [ "$RUN_NAME" = "GS1p1_tide_final" ]; then
  export FN_NEMO_DATA="/gws/nopw/j04/class_vol2/senemo/RUNS2024r01/GS1p1_tide_final/output/SENEMO_1y_19810101_19811231_grid_T_2D.nc"
  export FN_NEMO_DOMAIN="/gws/nopw/j04/class_vol2/senemo/RUNS2024r01/GS1p1_tide_final/config/domain_cfg.nc"

elif [ "$RUN_NAME" = "GS1p7_triads" ]; then
  export FN_NEMO_DATA="/gws/nopw/j04/class_vol2/senemo/jdha/GS1p7_TRIADS/output/SENEMO_1y_19810101_19811231_grid_T_2D.nc"
  export FN_NEMO_DOMAIN="/gws/nopw/j04/class_vol2/senemo/jdha/GS1p7_TRIADS/config/domain_cfg.nc"
fi

export GRID_OBS_RAD=30 # limit of acceptable distance between wet point (grid) and observation (km)
export THIN_OBS_RAD=50 # radius for thinning the observations (km)

#export FN_HARM_OBS="/home/users/jelt/data/obs/for_DA_dense/obs_M2.nc"  # Original obs data from GTM work
#export FN_HARM_OBS="/gws/nopw/j04/class_vol2/senemo/jelt/data/obs/for_validation_sparse/obs_M2.nc"  # Generated with submit_pre_process_obs.sh
export FN_HARM_OBS="/gws/nopw/j04/class_vol2/senemo/jelt/data/for_validation_"$GRID_OBS_RAD"_"$THIN_OBS_RAD"/obs_M2.nc"  # Generated with submit_pre_process_obs.sh


export FN_NEMO_CFG="/home/users/jelt/GitHub/COAsT/config/example_nemo_grid_t.json"

export DN_OUT="/gws/nopw/j04/class_vol2/senemo/jelt/"

export FN_ANALYSIS_OUT=$DN_OUT"PROCESSED_"$GRID_OBS_RAD"_"$THIN_OBS_RAD"/"$RUN_NAME"_extracted.nc"



