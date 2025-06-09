#export MACHINE="LOTUS"  # resource on JASMIN. Already set

#export CONDA_ENV="/home/users/jelt/miniconda3/envs/coast_dev"
export CONDA_ENV="/home/users/jelt/.conda/envs/coast_cmo_dev"
# location of COAsT repo, if using a particular branch
export COAST_REPO="/home/users/jelt/GitHub/COAsT"


#export FN_NEMO_DATA="/Users/jelt/Downloads/SENEMO/TIDE/SENEMO_1y_19810101_19811231_grid_T_2D.nc"
#export FN_NEMO_DOMAIN="/Users/jelt/Downloads/SENEMO/TIDE/domain_cfg.nc"
#export DN_FES="/Users/jelt/DATA/FES2014/ocean_tide_extrapolated/"
#export FN_OUT_DIR="/Users/jelt/Downloads/SENEMO/data/for_validation_sparse/obs_"
#export DN_OBS="/Users/jelt/GitHub/NEMO_validation/TG_obs_preprocessing/data/obs/"

export FN_NEMO_DATA="/gws/nopw/j04/class_vol2/senemo/jdha/GS1p7_TRIADS/output/SENEMO_1y_19810101_19811231_grid_T_2D.nc"
export FN_NEMO_DOMAIN="/gws/nopw/j04/class_vol2/senemo/jdha/GS1p7_TRIADS/config/domain_cfg.nc"
export DN_FES="/gws/nopw/j04/class_vol2/senemo/shared/FES2014/"
export FN_OUT_DIR="/gws/nopw/j04/class_vol2/senemo/jelt/data/for_validation_sparse/obs_"
export DN_OBS="/gws/nopw/j04/class_vol2/senemo/jelt/data/obs/"




#export FN_OUT_DIR="/home/users/jelt/data/obs/for_DA_dense/obs_M2.nc"  # Original obs data from GTM work
#export FN_OUT_DIR="/gws/nopw/j04/class_vol2/senemo/jelt/data/for_validation_sparse/obs_"  # Generated with pre_process_obs.py


export DN_OUT="/gws/nopw/j04/class_vol2/senemo/jelt/"

export FN_ANALYSIS_OUT=$DN_OUT"PROCESSED/"$RUN_NAME"_extracted.nc"



