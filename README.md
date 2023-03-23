# NEMO_validation
Scripts for validation of Coastal Ocean NEMO output.

This is split into a number of sections with corresponding directories:
```
EN4_preprocessing
EN4_processing
EN4_postprocessing
```

There are also original files that are under review (in `review_in_progress_Orig_files`, which includes sea level
analysis) and two directories that are temporary.

# Steps for processing data 

## Preprocess the EN4 data

Prior to computing diagnostics, the EN4 data is preprocessed into monthly COAsT-ready files for a restricted geographic domain.
To make it easier to work across sites and architectures two config files are used to control python and bash processes.

1. `cd EN4_preprocessing`

2. `config.sh` and `<MACHINE>_config.sh` must both be edited for machine choices, conda environment, paths etc.

Then preprocesing is triggered with:

3. Execute `. ./iter_en4_proc.sh`

which calls a machine dependant scheduler script `$MACHINE_pre_process_en4_monthly.sh` to invoke `pre_process_en4_monthly.py`

Output files are stored in the directory `config.sh: DOUT_EN4` with file structure: `<region>_processed_yyyymm.nc`. E.g.
`AMM15_processed_201410.nc`

## Process monthly data to generate model error profiles with respect to EN4 profiles

1. `cd EN4_processing`

2. `config.sh` and `<MACHINE>_config.sh` must both be edited for machine choices, conda environment, paths etc.

We use `iter_sub_METEST.sh`  to submit over all years and months separately. This allows for simple parallelisation 
as each month can be independently processed. This script sets the paths and variable names and launches a machine specific
script to process each month.

```
sbatch ${MACHINE,,}_ana_MOD_METEST.sh $MOD $start $month $end $GRID
```

where:

* $MOD is the Experiment e.g. P0.0
* $start is the start year
* $month is the month
* $end is the endyear
* $GRID contains is the domain file with grid info for that experiment

`spice_ana_MOD_METEST.sh` in turn calls the machine independent python script:

```
python  GEN_MOD_Dave_example_profile_validation.py $1 $2 $3 $4 $5  > LOGS/OUT_$1_$2_$3_$4_$5.log
```
using arguments: $1 $2 $3 $4 $5 corresponding to the above.

This outputs, in `DN_OUT/$REGION/`, files like: 
```
extracted_profiles_p0_200401_2005.nc
interpolated_profiles_p0_200401_2005.nc
interpolated_obs_p0_200401_2005.nc
profile_errors_p0_200401_2005.nc
surface_data_p0_200401_2005.nc
mid_data_p0_200401_2005.nc
bottom_data_p0_200401_2005.nc
mask_means_daily_p0_200401_2005.nc

```

However, some months have many profiles and some months few, so they take differing times to complete on different nodes.
Experience found that most months were completed in 20mins, about 10% needed 1hr, 5% 2hr and a couple needed 3hrs.
A short script with commandline control of the allocated walltime can see the slowest jobs, which previously ran out of 
walltime, through. For example:
```
#!/bin/bash
# comment out --time in lotus_ana_MOD_METEST.sh so it can be specified here
echo "Bash version ${BASH_VERSION}..."
source config.sh

rm LOGS/OUT* LOGS/*.err LOGS/*.out

#sbatch -J 201407 --time=2:00:00 lotus_ana_MOD_METEST.sh P0.0 2014 7 2015 CO7_EXACT_CFG_FILE.nc
#sbatch -J 201010 --time=2:00:00 lotus_ana_MOD_METEST.sh P0.0 2010 10 2011 CO7_EXACT_CFG_FILE.nc
#sbatch -J 201011 --time=2:00:00 lotus_ana_MOD_METEST.sh P0.0 2010 11 2011 CO7_EXACT_CFG_FILE.nc
sbatch -J 201109 --time=3:00:00 lotus_ana_MOD_METEST.sh P0.0 2011 9 2012 CO7_EXACT_CFG_FILE.nc
#sbatch -J 201110 --time=2:00:00 lotus_ana_MOD_METEST.sh P0.0 2011 10 2012 CO7_EXACT_CFG_FILE.nc
sbatch -J 200905 --time=3:00:00 lotus_ana_MOD_METEST.sh P0.0 2009 5 2010 CO7_EXACT_CFG_FILE.nc
```

## Postprocessing

1. `cd EN4_postprocessing`

2. `config.sh` and `<MACHINE>_config.sh` must both be edited for machine choices, conda environment, paths etc.

### Concatenate error profiles (merge seasons)

Merge seasons (DJF, MAM, JJA, SON) from multiple years into single files.

Execute with:
```
iter_merge_season.sh
```
which is just a simple concatenating loop over each season.

Each month invokes a machine specific sbatch scripts (e.g `spice_merge_season.sbatch`) where the model and season are 
passed onto a generic script
`python merge_season.py $1 $2 #1=Model, 2=month` 

Outputs are written to DN_OUT by season string, sss:
```
sss_PRO_INDEX.nc  ## merging interpolated_profiles_*.nc (model profiles on ref levels)
sss_PRO_DIFF.nc   ## merging profile_errors_*.nc (diff between model & obs on ref levels)
```

### Create Means

Then call `iter_mean_season.sh` to compute the spatial means over subregions within the NWS domain.

This launches machine specific script 

`sbatch ${MACHINE,,}_mean_season.sbatch $MOD $month`
that in turn launches a machine independent script:
```
python mean_season.py $1 $2 > LOGS/mean_season_$1_$2.log  # 1=Model, 2=month
```

to compute averages in each of the defined regions:
```
region_names = [ 'N. North Sea','S. North Sea','Eng. Channel','Outer Shelf', 'Irish Sea', 
                    'Kattegat', 'Nor. Trench', 'FSC', 'Off-shelf']
```
Creating output:
```
DJF_mask_means_daily.nc
MAM_mask_means_daily.nc
JJA_mask_means_daily.nc
SON_mask_means_daily.nc
```


### Plot the results.

Plot panels of regionally averaged profiles for summer and winter.

```
iter_plot_season.sh
```

sets machine specific paths and variables and launches

```sbatch ${MACHINE,,}_plot_season.sbatch```

which submits the following machine independent script

```python plot_season.py > LOGS/plot_season.log```

This plots multiple panels of area meaned profiles. One panel per region. Top row DJF and lower row JJA.
This also iterates over variables (temperature, salinity) and diagnostics (MAE, Bias).


Outputs e.g. `FIGS/regional_means_abs_diff_salinity_test.svg`
