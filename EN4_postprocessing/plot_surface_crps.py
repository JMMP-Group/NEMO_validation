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
analysis_str = "CRPS"

# List of analysis output files. Profiles from each will be plotted
# on each axis of the plot

# Assuming loading two configs: co7 and the P0.0. HARD WIRING. NOT IDEAL
fn_list = [config.dn_out+"All_mask_means_crps_co7.nc", config.dn_out.replace('co7','co7')+"All_mask_means_crps_co7.nc"]#

#%% General Plot Settings
## Specify the specific regions, their labels and order of appearance to plot. Region indexing to match EN4_postprocessing mean_season.py
region_ind = [ 1, 7, 3, 2, 9, 5, 4, 6, 8]              # Region indices (in analysis) to plot
region_names = [ 'N. North Sea','S. North Sea','Eng. Channel','Outer Shelf', 'Irish Sea', 'Kattegat', 'Nor. Trench', 'FSC', 'Off-shelf']

save_plot = True                # Boolean to save plot or not

# Subplot axes settings
n_r = 1               # Number of subplot rows
n_c = 9 #7               # Number of subplot columns
figsize = (9,4)       # Figure size
sharey = True         # Align y axes
sharex = False        # Align x axes
subplot_padding = .5  # Amount of vertical and horizontal padding between plots
fig_pad = (.075, .125, .12, .15)   # Whole figure padding as % (left, bottom, right, top)

max_crps_sal = 2.25
max_crps_temp = 1.2

# Legend settings, List of strings to use in legend (match with fn_list ordering)
legend_str = ["CO7","CO9p0"]
legend_index = 8 #11          # Axis index to put legend (flattened index, start from 0).
                          # Good to place in an empty subplot
legend_pos = 'upper right' # Position for legend, using matplitlib legend string
legend_fontsize = 10

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
 if var_str == "Temperature":
     max_crps = max_crps_temp
     units = "(deg C)"
 if var_str == "Salinity":
     max_crps = max_crps_sal
     units = "(psu)"


 if(1):
    # Patch in MAE for radius = 0
    mae_var_name = "profile_mean_abs_diff_{0}".format(var_str.lower())  # profile_mean_abs_diff_temperature



    # Labels and Titles
    xlabel = "{0} (units)".format(analysis_str)           # Xlabel string
    fig_title = "Regional {0} || {1}".format(analysis_str, var_str)  # Whole figure title


    var_name = "profile_mean_{0}_crps".format(var_str.lower())  # profile_mean_temperature_crps

    # Filename for the output
    fn_out = "FIGS/regional_{0}.svg".format(var_name)

    # Create plot and flatten axis array
    f,a = plt.subplots(n_r, n_c, figsize = figsize, sharey = sharey)
    a_flat = a.flatten()
    a_flat = a

    # Loop over regions
    for ii in range(n_reg):

        # Get the index of this region
        index = region_ind[ii]

        # Loop over datasets and plot their variable
        p = []
        for pp in range(n_ds):
            ds = ds_list[pp]
            # replace radius=0 with mae
            mae = ds[mae_var_name]
            rad, _ = xr.broadcast(ds.radius, ds[var_name])
            ds[var_name] = xr.where(rad == 0, mae.isel(dim_mask=index), ds[var_name])
            p.append( a_flat[ii].plot(ds.radius, ds[var_name].isel(dim_mask=index), '-o')[0])

        # Do some plot things
        a_flat[ii].set_title(f"{region_names[ii]}", fontsize=8)
        a_flat[ii].grid()
        a_flat[ii].set_ylim(0, max_crps)
        a_flat[ii].set_xlim(0, 20)

 if(1):
  # Make legend
  #a_flat[legend_index].legend(p, legend_str, loc=7, fontsize = legend_fontsize)
  f.legend(p, legend_str, fontsize = legend_fontsize,  bbox_to_anchor=(1.0,0.58))

  # Set Figure title and labels axes
  f.suptitle(fig_title, fontsize = title_fontsize, fontweight = title_fontweight)
  f.text(0.5, 0.01, 'radius (km)', ha='center')
  f.text(0.01, 0.5, 'CRPS' + " " + units, rotation=90, va='center')


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
