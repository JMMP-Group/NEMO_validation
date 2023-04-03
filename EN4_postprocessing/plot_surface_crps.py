'''
For plotting analysis data from a netcdf file created using COAsT.ProfileAnalysis.mask_means().
This will plot multiple datasets onto a set of subplots. Each subplot is for
a different averaging region.

At the top of this script, you can set the paths to the netcdf files to plot
and where to save. If you have multiple model runs to plot, provide a list
of file paths (strings).

Below this section are a bunch of parameters you can set, with explanations in
comments. Edit this as much as you like or even go into the plotting code below.
'''

from config import config
config = config() # initialise variables in python

config.dn_out = "/Users/jelt/Downloads/"

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

#%% File settings
run_name = "test"


# List of analysis output files. Profiles from each will be plotted
# on each axis of the plot

# Assuming loading two configs: co7 and the P0.0. HARD WIRING. NOT IDEAL
fn_list = [config.dn_out+"All_mask_means_crps_co7.nc", config.dn_out.replace('co7','co7')+"All_mask_means_crps_co7.nc"]#

#%% General Plot Settings
## Specify the specific regions, their labels and order of appearance to plot. Region indexing to match EN4_postprocessing mean_season.py
region_ind = [ 1, 7, 3, 2, 9, 5, 4, 6, 8]              # Region indices (in analysis) to plot
region_names = [ 'N. North Sea','S. North Sea','Eng. Channel','Outer Shelf', 'Irish Sea', 'Kattegat', 'Nor. Trench', 'FSC', 'Off-shelf']

plot_zero_line = False            # Plot a black vertical line at x = 0
plot_mean_depth = False          # Plot the mean bathymetric depth. Make sure 'bathymetry' is in the analysis dataset
save_plot = True                # Boolean to save plot or not

#ref_depth = np.concatenate((np.arange(1,100,2), np.arange(100,300,5), np.arange(300, 1000, 50))) # Data depths
#ref_depth = np.concatenate((np.arange(1,100,2), np.arange(100,300,5), np.arange(300, 1000, 50), np.arange(1000,2000,100)))
# Should match definition in EN4_processing: ref_depth
ref_depth = np.concatenate((np.arange(1,100,2), np.arange(100,300,5), np.arange(300, 1000, 50), np.arange(1000,4000,100)))
#ref_depth = np.arange(1,4000,10)

print (np.shape(ref_depth))

# Subplot axes settings
n_r = 1               # Number of subplot rows
n_c = 9 #7               # Number of subplot columns
figsize = (9,7)       # Figure size
sharey = True         # Align y axes
sharex = False        # Align x axes
subplot_padding = .5  # Amount of vertical and horizontal padding between plots
fig_pad = (.075, .075, .075, .1)   # Whole figure padding as % (left, top, right, bottom)
#max_depth = 50                    # Maximum plot depth
#max_depth = 4000                   # Maximum plot depth
#max_depth = 800                   # Maximum plot depth
max_depth = 150                   # Maximum plot depth
max_crps_sal = 2.25
max_crps_temp = 1.0

# Legend settings
#legend_str = ["P0.0", "P0.9", "BOB"]     # List of strings to use in legend (match with fn_list ordering)
#legend_str = ["P0.0","P0.1b","P0.6","u-cp812","u-cq846" ]     # List of strings to use in legend (match with fn_list ordering)
#legend_str = ["P0.0","GLS","ERA5","GLS+ERA5","GEG_TPX_SF12","GEG_FES_SF12","GEG_FES_ME" ]     # List of strings to use in legend (match with fn_list ordering)
#legend_str = ["P0.0","GLS","ERA5","GLS+ERA5","P0.0NERCWAD","P0.0NERCNOWAD","GEG_FES_ME"]# ,"GEG_TPX_SF12","GEG_FES_SF12","GEG_FES_ME" ]     # List of strings to use in legend (match with fn_list ordering)
#legend_str = ["P0.0","P0.0NERCWAD","P0.0NERCNOWAD","p0 GEG OLDBDY 0.69", "p0 GEG OLD BDY 0.7","P0.0_NERC_GEG"]# ,"GEG_TPX_SF12","GEG_FES_SF12","GEG_FES_ME" ]     # List of strings to use in legend (match with fn_list ordering)
#legend_str = ["P0.0","P0.6","P0.7","P0.8","503","504","535","545","P0.5.c"]# ,"GEG_TPX_SF12","GEG_FES_SF12","GEG_FES_ME" ]     # List of strings to use in legend (match with fn_list ordering)
legend_str = ["CO7"] #,"P0.0"]
#legend_str = ["JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","OCT","NOV","DEC"]# ,"GEG_TPX_SF12","GEG_FES_SF12","GEG_FES_ME" ]     # List of strings to use in legend (match with fn_list ordering)
#legend_str = ["P0.0","P0.1b","P0.6","P0.7","P0.8","u-cp812","u-cp815","u-cq846" ]     # List of strings to use in legend (match with fn_list ordering)
#legend_str = ["P0.0","FES_P0.1b","TPX_GLS_P0.6","TPX_ERA5_P0.7","TPX_GLS_ERA5_P0.8","FES_MEu-cq846" ]     # List of strings to use in legend (match with fn_list ordering)
legend_index = 8 #11          # Axis index to put legend (flattened index, start from 0).
                          # Good to place in an empty subplot
legend_pos = 'upper right' # Position for legend, using matplitlib legend string
legend_fontsize =  6

# Labels and Titles
xlabelpos = (figsize[0]/2, 0)              # (x,y) position of xlabel
ylabel = "Depth (m)"                       # Ylabel string
ylabelpos = (figsize[1]/2, 0)              # (x,y) position of ylabel

label_fontsize = 8 #11                        # Fontsize of all labels
label_fontweight = "normal"                # Fontweight to use for labels and subtitles
title_fontsize = 13                        # Fontsize of title
title_fontweight = "bold"                  # Fontweight to use for title




#%% SCRIPT: READ AND PLOT DATA

# Read all datasets into list
ds_list = [xr.open_dataset(dd) for dd in fn_list]
n_ds = len(ds_list)
n_reg = len(region_ind)
n_ax = n_r*n_c

print(f"Check region names specified are consistent with mask file")
for i in range(n_reg):
  print(f"Panel label:({region_names[i]}) matches data label: ({ds_list[0].region_names.values[region_ind[i]]})")

# Loop over variable to plot
for var_str in ["Temperature", "Salinity"]:
 if var_str == "Temperature": max_crps = max_crps_temp
 if var_str == "Salinity": max_crps = max_crps_sal

 for analysis_str in ["CRPS"]:

  # Patch in MAE for radius = 0
  mae_var_name = "profile_mean_abs_diff_{0}".format(var_str.lower())  # profile_mean_abs_diff_temperature



  # Labels and Titles
  xlabel = "{0} (units)".format(analysis_str)           # Xlabel string
  fig_title = "Regional {0} || {1}".format(analysis_str, var_str)  # Whole figure title


  var_name = "profile_mean_{0}_crps".format(var_str.lower())  # profile_mean_temperature_crps

  # Filename for the output
  fn_out = "FIGS/regional_{0}_{1}.svg".format(var_name, run_name)

  # Create plot and flatten axis array
  f,a = plt.subplots(n_r, n_c, figsize = figsize, sharex = sharex, sharey = sharey)
  a_flat = a.flatten()
  a_flat = a

  # Loop over regions
  #for ii in range(n_ax):
  for ii in range(n_reg):
    print (ii,n_ax)
    for row, season in enumerate(["DJF"]):
        if ii >= n_reg:
            a_flat[row,ii].axis('off')
            continue

        # Get the index of this region
        index = region_ind[ii]
	
        # Loop over datasets and plot their variable
        p = []
        for pp in range(n_ds):
            ds = ds_list[pp]
            mae = ds[mae_var_name]
            rad, _ = xr.broadcast(ds.radius, ds[var_name])
            ds[var_name] = xr.where(rad == 0, mae.isel(dim_mask=index), ds[var_name])
            print(ds[var_name])

            p.append( a_flat[ii].plot(ds.radius, ds[var_name].isel(dim_mask=index), c='b')[0] )
            p.append( a_flat[ii].plot(ds.radius, ds[var_name].isel(dim_mask=index), 'bo')[0])

        # Do some plot things
        a_flat[ii].set_title(f"{region_names[ii]}", fontsize=8)
        a_flat[ii].grid()
        a_flat[ii].set_ylim(0, max_crps)
        a_flat[ii].set_xlim(0, 20)

  # Make legend
  a_flat[legend_index].legend(p, legend_str, loc=7, fontsize = legend_fontsize)

  # Set Figure title and labels axes
  f.suptitle(fig_title, fontsize = title_fontsize, fontweight = title_fontweight)
  f.text(0.5, 0.01, 'radius (km)')
  f.text(0.01, 0.5, 'CRPS', rotation=90)


  # Set x and y labels
  f.text(xlabelpos[0], xlabelpos[1], xlabel, 
       va = 'center', rotation = 'horizontal', 
       fontweight = label_fontweight, fontsize = label_fontsize)

  # Set tight figure layout
  f.tight_layout(w_pad = subplot_padding, h_pad= subplot_padding)
  f.subplots_adjust(left   = (fig_pad[0]), 
                  bottom = (fig_pad[1]),
                  right  = (1-fig_pad[2]),
                  top    = (1-fig_pad[3]))

  # Save plot maybe
  if save_plot: 
      f.savefig(fn_out)
