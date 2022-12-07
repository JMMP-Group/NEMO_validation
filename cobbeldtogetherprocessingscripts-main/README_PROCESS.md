# Steps for procesing data 

# Process the en4 data

before we do anythgin we break up the en4 data into processed months for coast

we use

* iter_en4_proc.sh to loop through and call
* pre_process_en4_monthly.py

## Process monthly data
In recent coast development env (my conda env coast_nov2022)



we use **iter_sub_METEST.sh**  to submit over all years and months separately,
its a bit dumb as some months have only a few profiles asn some much more so the queue time 
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


## Concatenate profiles (merge months)

**The next steps I used my older version of COAST, but should be easy enough to convert**

**iter_merge.sh** is just a simple loop over months

calling **merge.sbatch** which in turn calls 

```bash
python merge_monthly.py $1 $2 #(Model month)**
```

That just concatenates all matching a specified month into a single file
in preparation for creating a mean for each month.


## Create Means

As the runs were not done yet  I just did them by hand on the command line but could have a  simple loop 
through the months as above

```bash
sbatch mean.sbatch P0.0 2
```
for example will mean up all for feb


## Plot the results.

To plot the months I used:
**plot_month.py**

This reads in each meaned month and plots the various regions.

To cut off say top data for the salinity I can contrain by index:

```python
p.append( a_flat[ii].plot(ds[var_name][index][4:150], ref_depth[4:150])[0] )
```

where the above starts at 4 instead of 0, indeed it also chops off the bottom levels only going to 150


 


