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

**The next steps I used my older version of COAST, but should be easy enough to convert**

**iter_merge.sh** is just a simple loop over months

calling **merge.sbatch** which in turn calls 

```bash
python merge_monthly.py $1 $2 #(Model month)**
```

That just concatenates all matching a specified month into a single file
in preparation for creating a mean for each month.


### Create Means

As the runs were not done yet  I just did them by hand on the command line but could have a  simple loop 
through the months as above

```bash
sbatch spice_mean.sbatch P0.0 2
```
for example will mean up all for feb


### Plot the results.

To plot the months I used:
**plot_month.py**

This reads in each meaned month and plots the various regions.

To cut off say top data for the salinity I can contrain by index:

```python
p.append( a_flat[ii].plot(ds[var_name][index][4:150], ref_depth[4:150])[0] )
```

where the above starts at 4 instead of 0, indeed it also chops off the bottom levels only going to 150
