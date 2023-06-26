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

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

#%% File settings
run_name = "test"

# List of analysis output files. Profiles from each will be plotted
# on each axis of the plot

# Assuming loading two configs: co7 and the P0.0. HARD WIRING. NOT IDEAL
co7_path = '/gws/nopw/j04/jmmp/CO9_AMM15_validation/co7/profiles/'
fn_list_DJF = [config.dn_out+"DJF_mask_means_daily.nc",
               co7_path+"DJF_mask_means_daily.nc"]#
fn_list_JJA = [config.dn_out+"JJA_mask_means_daily.nc",
               co7_path+"JJA_mask_means_daily.nc"]#
print (fn_list_JJA)
#fn_list_DJF = [config.dn_out+"DJF_mask_means_daily.nc"] 
#fn_list_JJA = [config.dn_out+"JJA_mask_means_daily.nc"]

#%% General Plot Settings
## Specify the specific regions, their labels and order of appearance to plot. Region indexing to match EN4_postprocessing mean_season.py
region_ind = [ 1, 7, 3, 2, 9, 5, 4, 6, 8]              # Region indices (in analysis) to plot
region_names = [ 'N. North Sea','S. North Sea','Eng. Channel','Outer Shelf', 'Irish Sea', 'Kattegat', 'Nor. Trench', 'FSC', 'Off-shelf']

#region_names = ["A","B","C","D","E","F","G","H","I"]  # Region names, will be used for titles in plot

plot_zero_line = True            # Plot a black vertical line at x = 0
plot_mean_depth = True          # Plot the mean bathymetric depth. Make sure 'bathymetry' is in the analysis dataset
save_plot = True                # Boolean to save plot or not

#ref_depth = np.concatenate((np.arange(1,100,2), np.arange(100,300,5), np.arange(300, 1000, 50))) # Data depths
#ref_depth = np.concatenate((np.arange(1,100,2), np.arange(100,300,5), np.arange(300, 1000, 50), np.arange(1000,2000,100)))
# Should match definition in EN4_processing: ref_depth
ref_depth = np.concatenate((np.arange(1,100,2), np.arange(100,300,5), np.arange(300, 1000, 50), np.arange(1000,4000,100)))
#ref_depth = np.arange(1,4000,10)

print (np.shape(ref_depth))

# Subplot axes settings
n_r = 2               # Number of subplot rows
n_c = 9 #7               # Number of subplot columns
figsize = (9,7)       # Figure size
sharey = True         # Align y axes
sharex = False        # Align x axes
subplot_padding = .5  # Amount of vertical and horizontal padding between plots
fig_pad = (.075, .075, .15, .1)   # Whole figure padding as % (left, top, right, bottom)
#max_depth = 50                    # Maximum plot depth
#max_depth = 4000                   # Maximum plot depth
#max_depth = 800                   # Maximum plot depth
max_depth = 150                   # Maximum plot depth

# Legend settings
#legend_str = ["P0.0", "P0.9", "BOB"]     # List of strings to use in legend (match with fn_list ordering)
#legend_str = ["P0.0","P0.1b","P0.6","u-cp812","u-cq846" ]     # List of strings to use in legend (match with fn_list ordering)
#legend_str = ["P0.0","GLS","ERA5","GLS+ERA5","GEG_TPX_SF12","GEG_FES_SF12","GEG_FES_ME" ]     # List of strings to use in legend (match with fn_list ordering)
#legend_str = ["P0.0","GLS","ERA5","GLS+ERA5","P0.0NERCWAD","P0.0NERCNOWAD","GEG_FES_ME"]# ,"GEG_TPX_SF12","GEG_FES_SF12","GEG_FES_ME" ]     # List of strings to use in legend (match with fn_list ordering)
#legend_str = ["P0.0","P0.0NERCWAD","P0.0NERCNOWAD","p0 GEG OLDBDY 0.69", "p0 GEG OLD BDY 0.7","P0.0_NERC_GEG"]# ,"GEG_TPX_SF12","GEG_FES_SF12","GEG_FES_ME" ]     # List of strings to use in legend (match with fn_list ordering)
#legend_str = ["P0.0","P0.6","P0.7","P0.8","503","504","535","545","P0.5.c"]# ,"GEG_TPX_SF12","GEG_FES_SF12","GEG_FES_ME" ]     # List of strings to use in legend (match with fn_list ordering)
legend_str = ["CO9p2","CO7"]
#legend_str = ["JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","OCT","NOV","DEC"]# ,"GEG_TPX_SF12","GEG_FES_SF12","GEG_FES_ME" ]     # List of strings to use in legend (match with fn_list ordering)
#legend_str = ["P0.0","P0.1b","P0.6","P0.7","P0.8","u-cp812","u-cp815","u-cq846" ]     # List of strings to use in legend (match with fn_list ordering)
#legend_str = ["P0.0","FES_P0.1b","TPX_GLS_P0.6","TPX_ERA5_P0.7","TPX_GLS_ERA5_P0.8","FES_MEu-cq846" ]     # List of strings to use in legend (match with fn_list ordering)
legend_index = 6 #11          # Axis index to put legend (flattened index, start from 0).
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
ds_list_DJF = [xr.open_dataset(dd) for dd in fn_list_DJF]
ds_list_JJA = [xr.open_dataset(dd) for dd in fn_list_JJA]
ds_list = ds_list_DJF
n_ds = len(ds_list)
n_reg = len(region_ind)
n_ax = n_r*n_c

print(f"Check region names specified are consistent with mask file")
for i in range(n_reg):
  print(f"Panel label:({region_names[i]}) matches data label: ({ds_list[0].region_names.values[region_ind[i]]})")

# Loop over variable to plot
for var_str in ["Temperature", "Salinity"]:
 if var_str == "Temperature": units = "(deg C)"
 if var_str == "Salinity": units = "(psu)"

 #for analysis_str in ["MAE", "STD", "BIAS"]:
 for analysis_str in ["MAE", "BIAS"]:
  # Labels and Titles
  xlabel = "{0} (units)".format(analysis_str)           # Xlabel string
  fig_title = "Regional {0} || {1}".format(analysis_str, var_str)  # Whole figure title
 
  if analysis_str == "MAE":
    tmp_str = "mean_abs_diff"
  if analysis_str == "STD":
    tmp_str = "std_diff"
  if analysis_str == "BIAS":
    tmp_str = "mean_diff"

  var_name = "profile_{0}_{1}".format(tmp_str, var_str.lower())

  # Filename for the output
  fn_out = "FIGS/regional_{0}_{1}.svg".format(var_name, run_name)

  # Create plot and flatten axis array
  f,a = plt.subplots(n_r, n_c, figsize = figsize, sharex = sharex, sharey = sharey)
  a_flat = a

  # Loop over regions
  #for ii in range(n_ax):
  for ii in range(n_reg):
    print (ii,n_ax)
    for row, season in enumerate(["DJF", "JJA"]):
        if ii >= n_reg:
            a_flat[row,ii].axis('off')
            continue

        # Get the index of this region
        index = region_ind[ii]
	
        # Loop over datasets and plot their variable
        p = []
        for pp in range(n_ds):

            print(f"season:{season}")
            if season == "DJF":
              ds = ds_list_DJF[pp]
            elif season == "JJA":
              ds = ds_list_JJA[pp]
            else:
              print(f"Not expecting that season: {season}")
            #p.append( a_flat[ii].plot(ds[var_name][index], ref_depth)[0] )
            #p.append( a_flat[ii].plot(ds[var_name][index][4:50], ref_depth[4:50])[0] )
            #p.append( a_flat[row,ii].plot(ds[var_name][index][4:150], ref_depth[4:150])[0] )
            #p.append( a_flat[row,ii].plot(ds[var_name][index][1:150], ref_depth[1:150])[0] )
            p.append( a_flat[row,ii].plot(ds[var_name][index][:100], ref_depth[:100])[0] )

        # Do some plot things
        a_flat[row,ii].set_title(f"{region_names[ii]}:\n{season}", fontsize=8)
        a_flat[row,ii].grid()
        a_flat[row,ii].set_ylim(0, max_depth)

        # set x lims
        if var_str == 'Salinity':
            a_flat[row,ii].set_xlim(-0.1, 0.5)
        if var_str == 'Temperature':
            a_flat[row,ii].set_xlim(-0.1, 2.5)

        # Plot fixed lines at 0 and mean depth
        if plot_zero_line:
            a_flat[row,ii].plot([0,0], [0, max_depth], c='k', linewidth = 1, linestyle = '-')
        if plot_mean_depth:
            a_flat[row,ii].plot([a_flat[row,ii].get_xlim()[0], a_flat[row,ii].get_xlim()[1]], 
                                [ds['bathymetry'][index],ds['bathymetry'][index]],
                                color='k', linestyle='--')

        # Invert y axis
        a_flat[row,ii].invert_yaxis()

  # Make legend
  #a_flat[0,legend_index].legend(p, legend_str, fontsize = legend_fontsize, bbox_to_anchor=(1.2, 0.6))
  f.legend(p, legend_str, fontsize = legend_fontsize,  bbox_to_anchor=(1.0,0.5))


  # Set Figure title and label axes
  f.suptitle(fig_title, fontsize = title_fontsize, fontweight = title_fontweight)
  f.text(0.5, 0.01, var_str + " " + units, ha='center')
  f.text(0.01, 0.5, 'depth (m)', va='center', rotation=90)

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
