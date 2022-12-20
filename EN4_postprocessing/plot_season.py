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
#fn_list = ["/Users/dbyrne/transfer/mask_means_daily_test.nc",
#           "/Users/dbyrne/transfer/mask_means_daily_test.nc"]
#fn_list = ["/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/mi-bd207/analysis/mask_means_daily_p0_2003_2004.nc", "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/rosie_mi-an561_1990/analysis/mask_means_daily_p0_2003_2004.nc"]
#fn_list = ["/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/P0.0/analysis/ALL_mask_means_daily.nc","/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/P0.9/analysis/ALL_mask_means_daily.nc","/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/P0.0/analysis/ALL_mask_means_daily.nc"]
fn_list = ["%s%03s_mask_means_daily.nc"%(config.dn_out, "DJF"),
           "%s%03s_mask_means_daily.nc"%(config.dn_out, "JJA")]#
tt=[           "%s%02d_mask_means_daily.nc"%(config.dn_out, 3),
           "%s%02d_mask_means_daily.nc"%(config.dn_out, 4),
           "%s%02d_mask_means_daily.nc"%(config.dn_out, 5),
           "%s%02d_mask_means_daily.nc"%(config.dn_out, 6),
           "%s%02d_mask_means_daily.nc"%(config.dn_out, 7),
           "%s%02d_mask_means_daily.nc"%(config.dn_out, 8),
           "%s%02d_mask_means_daily.nc"%(config.dn_out, 9),
           "%s%02d_mask_means_daily.nc"%(config.dn_out, 10),
           "%s%02d_mask_means_daily.nc"%(config.dn_out, 11),
           "%s%02d_mask_means_daily.nc"%(config.dn_out, 12)]
#           "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/P0.0/analysisb/01_mask_means_daily.nc",
#           "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/P0.0/analysisb/02_mask_means_daily.nc",
#           "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/P0.0/analysisb/03_mask_means_daily.nc",
#           "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/P0.0/analysisb/04_mask_means_daily.nc",
#           "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/P0.0/analysisb/06_mask_means_daily.nc",
#           "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/P0.0/analysisb/07_mask_means_daily.nc",
#           "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/P0.0/analysisb/08_mask_means_daily.nc",
#           "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/P0.0/analysisb/10_mask_means_daily.nc",
#           "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/P0.0/analysisb/11_mask_means_daily.nc",
#           "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/P0.0/analysisb/12_mask_means_daily.nc"]
#           "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/P0.6/analysis/ALL_mask_means_daily.nc",
#           "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/P0.7/analysis/ALL_mask_means_daily.nc",
#           "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/P0.8/analysis/ALL_mask_means_daily.nc",
#           "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/u-cr503/analysis/ALL_mask_means_daily.nc",
#           "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/u-cr504/analysis/ALL_mask_means_daily.nc",
#           "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/u-cr535/analysis/ALL_mask_means_daily.nc",
#           "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/u-cr545/analysis/ALL_mask_means_daily.nc",
#           "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/u-cr934/analysis/ALL_mask_means_daily.nc",
#           "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/P0.5.c/analysis/ALL_mask_means_daily.nc"]
#           "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/u-cq846/analysis/ALL_mask_means_daily.nc"]
#           "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/u-cp812/analysis/ALL_mask_means_daily.nc"]
##          "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/u-cp815/analysis/ALL_mask_means_daily.nc",
#          "/scratch/fred/COMPARE_VN36_VN_4.0_TIDE_SSH/u-cq846/analysis/ALL_mask_means_daily.nc"]#

# Filename for the output
fn_out = "FIGS/regional_means_{0}.svg".format(run_name)

#%% General Plot Settings
# regions need to match those in EN4_postprocessing mean_monthly.py
#region_ind = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]  # Region indices (in analysis) to plot
#region_names = ['whole_domain', 'N_north_sea', 'outer_shelf', 'eng_channel', 'nor_trench',
#                'kat', 'fsc', 'S_north_sea', 'off-shelf', 'Irish_Sea']
region_ind = [ 1, 7, 3, 2, 9, 5, 4]              # Region indices (in analysis) to plot
region_names = [ 'N. North Sea','S. North Sea','Eng. Channel','Outer Shelf', 'Irish Sea', 'Kattegat', 'Nor. Trench']

#region_names = ["A","B","C","D","E","F","G","H","I"]  # Region names, will be used for titles in plot
#var_name = "profile_mean_diff_temperature"     # Variable name in analysis file to plot
#var_name = "profile_mean_diff_salinity"     # Variable name in analysis file to plot
var_name = "profile_mean_abs_diff_temperature"     # Variable name in analysis file to plot
#var_name = "profile_mean_abs_diff_salinity"     # Variable name in analysis file to plot
plot_zero_line = True            # Plot a black vertical line at x = 0
plot_mean_depth = False          # Plot the mean bathymetric depth. Make sure 'bathymetry' is in the analysis dataset
save_plot = True                # Boolean to save plot or not

#ref_depth = np.concatenate((np.arange(1,100,2), np.arange(100,300,5), np.arange(300, 1000, 50))) # Data depths
#ref_depth = np.concatenate((np.arange(1,100,2), np.arange(100,300,5), np.arange(300, 1000, 50), np.arange(1000,2000,100)))
# Should match definition in EN4_processing: ref_depth
ref_depth = np.concatenate((np.arange(1,100,2), np.arange(100,300,5), np.arange(300, 1000, 50), np.arange(1000,4000,100)))
#ref_depth = np.arange(1,4000,10)

print (np.shape(ref_depth))

# Subplot axes settings
n_r = 2               # Number of subplot rows
n_c = 7               # Number of subplot columns
figsize = (7,7)       # Figure size
sharey = True         # Align y axes
sharex = False        # Align x axes
subplot_padding = .5  # Amount of vertical and horizontal padding between plots
fig_pad = (.075, .075, .075, .1)   # Whole figure padding as % (left, top, right, bottom)
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
legend_str = ["P0.0","P0.6","P0.7","P0.8","503","504","535","545","P0.5.c"]# ,"GEG_TPX_SF12","GEG_FES_SF12","GEG_FES_ME" ]     # List of strings to use in legend (match with fn_list ordering)
#legend_str = ["JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","OCT","NOV","DEC"]# ,"GEG_TPX_SF12","GEG_FES_SF12","GEG_FES_ME" ]     # List of strings to use in legend (match with fn_list ordering)
#legend_str = ["P0.0","P0.1b","P0.6","P0.7","P0.8","u-cp812","u-cp815","u-cq846" ]     # List of strings to use in legend (match with fn_list ordering)
#legend_str = ["P0.0","FES_P0.1b","TPX_GLS_P0.6","TPX_ERA5_P0.7","TPX_GLS_ERA5_P0.8","FES_MEu-cq846" ]     # List of strings to use in legend (match with fn_list ordering)
legend_index = 11          # Axis index to put legend (flattened index, start from 0).
                          # Good to place in an empty subplot
legend_pos = 'upper right' # Position for legend, using matplitlib legend string
legend_fontsize =  6

# Labels and Titles
xlabel = "Absolute Error (degC)"           # Xlabel string
xlabelpos = (figsize[0]/2, 0)              # (x,y) position of xlabel
ylabel = "Depth (m)"                       # Ylabel string
ylabelpos = (figsize[1]/2, 0)              # (x,y) position of ylabel
fig_title = "Regional MAE || All Seasons"  # Whole figure title
#fig_title = "Regional Mean Error "  # Whole figure title
label_fontsize = 11                        # Fontsize of all labels
label_fontweight = "normal"                # Fontweight to use for labels and subtitles
title_fontsize = 13                        # Fontsize of title
title_fontweight = "bold"                  # Fontweight to use for title




#%% SCRIPT: READ AND PLOT DATA

# Read all datasets into list
ds_list = [xr.open_dataset(dd) for dd in fn_list]
n_ds = len(ds_list)
n_reg = len(region_ind)
n_ax = n_r*n_c

# Create plot and flatten axis array
f,a = plt.subplots(n_r, n_c, figsize = figsize, sharex = sharex, sharey = sharey)
a_flat = a.flatten()

# Loop over regions
for ii in range(n_ax):
    print (ii,n_ax)
    
    if ii >= n_reg:
        a_flat[ii].axis('off')
        continue
    
    # Get the index of this region
    index = region_ind[ii]
    
    # Loop over datasets and plot their variable
    p = []
    for pp in range(n_ds):
        ds = ds_list[pp]
        print("ii is ",ii)
        print("size is" ,len(a_flat[:]) )
        print("size ds" ,len(ds[var_name][:]) )
        print("index is" ,index )
        #p.append( a_flat[ii].plot(ds[var_name][index], ref_depth)[0] )
        #p.append( a_flat[ii].plot(ds[var_name][index][4:50], ref_depth[4:50])[0] )
        #p.append( a_flat[ii].plot(ds[var_name][index][4:150], ref_depth[4:150])[0] )
        p.append( a_flat[ii].plot(ds[var_name][index][:100], ref_depth[:100])[0] )
        
    # Do some plot things
    a_flat[ii].set_title(region_names[ii])
    a_flat[ii].grid()
    a_flat[ii].set_ylim(0, max_depth)
        
    # Plot fixed lines at 0 and mean depth
    if plot_zero_line:
        a_flat[ii].plot([0,0], [0, max_depth], c='k', linewidth = 1, linestyle = '-')
    if plot_mean_depth:
        a_flat[ii].plot()
        
    # Invert y axis
    a_flat[ii].invert_yaxis()

# Make legend
a_flat[legend_index].legend(p, legend_str, fontsize = legend_fontsize)

# Set Figure title
f.suptitle(fig_title, fontsize = title_fontsize, fontweight = title_fontweight)

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
