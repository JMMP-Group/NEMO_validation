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

1 . `config.sh` must both be edited for path, conda environment, and machine choices.

Then preprocesing is triggered with:

2. Execute `. ./iter_en4_proc.sh`

which calls a machine dependant scheduler script `$MACHINE_pre_process_en4_monthly.sh` to invoke `pre_process_en4_monthly.py`

Output files are stored in the directory `config.sh: DOUT_EN4` with file structure: `<region>_processed_yyyymm.nc`

## Process monthly data to generate model error profiles with respect to EN4 profiles

We use **iter_sub_METEST.sh**  to submit over all years and months separately,
its a bit dumb as some months have only a few profiles and some much more so the queue time 
for the serial sbatch should really be a bit more dynamic e.g. the longer ones take hours
the shorter ones take minutes

anyway it calls:

**spice_ana_MOD_METEST.sh**

sending the arguments $MOD $start $month $end $GRID

where:

* $MOD is the Experiment e.g. P0.0
* $start is the start year
* $month is the month
* $end is the endyear
* $GRID contains is the domain file with grid info for that experiment

**spice_ana_MOD_METEST.sh** in turn calls the python 

**GEN_MOD_Dave_example_profile_validation.py**

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

## Postprocessing

### Concatenate error profiles (merge months)

Merge months from multiple years into single files.

Execute with:
```
iter_merge.sh
```
which is just a simple concatenating loop over each month.

Each month invokes a machine specific sbatch scripts (e.g `spice_merge.sbatch`) where the model and month are passed
onto a generic script
`python merge_monthly.py $1 $2 #1=Model, 2=month` 

Outputs are written to DN_OUT by month number, mm:
```
mm_PRO_INDEX.nc  ## merging interpolated_profiles_*.nc (model profiles on ref levels)
mm_PRO_DIFF.nc   ## merging profile_errors_*.nc (diff between model & obs on ref levels)
```

### Create Means

Then call `iter_mean.sh` to compute the spatial means over subregions within the NWS domain.

This launches machine specific script 

`sbatch ${MACHINE,,}_mean.sbatch $MOD $month`
that in turn launches a machine independent script:
```
python mean_monthly.py $1 $2 > LOGS/mean_monthly_$1_$2.log  # 1=Model, 2=month
```

to compute averages in each of the defined regions:
```
masks_names = ['whole_domain', 'north_sea','outer_shelf','eng_channel','nor_trench',
                    'kattegat','fsc','shelf_break', 'southern_north_sea', 'irish_sea' ]
```
Creating output:
```
01_mask_means_daily.nc
...
12_mask_means_daily.nc
```





### Plot the results.
#### JP:

```
iter_plot.sh
```

sets machine specific paths and variables and launches

```sbatch ${MACHINE,,}_plot.sbatch```

which submits the following machine independent script

```python plot_month.py > LOGS/plot_month.log```

This plots multiple panels of area meaned profiles. One panel per region.
The month, variable and depth properties are currently hardwired into `plot_month.py`


Outputs `FIGS/regional_means_test.svg`

#### WIP - Enda:
To plot the months I used:
**plot_month.py**

This reads in each meaned month and plots the various regions.

To cut off say top data for the salinity I can contrain by index:

```python
p.append( a_flat[ii].plot(ds[var_name][index][4:150], ref_depth[4:150])[0] )
```

where the above starts at 4 instead of 0, indeed it also chops off the bottom levels only going to 150
