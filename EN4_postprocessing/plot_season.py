from PythonEnvCfg.config import config
config = config() # initialise variables in python

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

class seasonal_profiles(object):
    '''
    For plotting analysis data from a netcdf file created using
    COAsT.ProfileAnalysis.mask_means(). 
    This will plot multiple datasets onto a set of subplots.
    Each subplot is for a different averaging region.
    
    At the top of this script, you can set the paths to the netcdf files to
    plot and where to save. If you have multiple model runs to plot, provide
    a list of file paths (strings).
    
    Below this section are a bunch of parameters you can set, with
    explanations in comments. 
    Edit this as much as you like or even go into the plotting code below.
    '''

    def __init__(self):
        #%% File settings
        self.run_name = "test"
        
        # List of analysis output files. Profiles from each will be plotted
        # on each axis of the plot
        
        # Assuming loading two configs: co7 and the P0.0. HARD WIRING.
        # NOT IDEAL
        co7_path = '/gws/nopw/j04/jmmp/CO9_AMM15_validation/co7/profiles/'
        self.fn_list_DJF = [config.dn_out+"profiles/DJF_mask_means_daily.nc",
                            co7_path+"DJF_mask_means_daily.nc"]
        self.fn_list_MAM = [config.dn_out+"profiles/MAM_mask_means_daily.nc",
                            co7_path+"MAM_mask_means_daily.nc"]
        self.fn_list_JJA = [config.dn_out+"profiles/JJA_mask_means_daily.nc",
                            co7_path+"JJA_mask_means_daily.nc"]
        self.fn_list_SON = [config.dn_out+"profiles/SON_mask_means_daily.nc",
                            co7_path+"SON_mask_means_daily.nc"]
    
    def plot_all_djf_jja(self):
        """
        Plot djf and jja for all locations.
        """

        #%% General Plot Settings
        ## Specify the specific regions, their labels and order of appearance
        ## to plot.
        ## Region indexing to match EN4_postprocessing mean_season.py

        # Region indices (in analysis) to plot
        self.region_ind = [ 1, 7, 3, 2, 9, 5, 4, 6, 8]
        self.region_names = ['N. North\nSea','S. North\nSea',
                        'Eng.\nChannel','Outer\nShelf',
                        'Irish\nSea', 'Kattegat',
                        'Nor.\nTrench', 'FSC', 'Off-shelf']
        
        self.plot_zero_line  = True      # Plot a black vertical line at x = 0
        # Plot the mean bathymetric depth. 
        # Make sure 'bathymetry' is in the analysis dataset
        self.plot_mean_depth = True      
        self.save_plot       = True      # Boolean to save plot or not
        
        # Should match definition in EN4_processing: ref_depth
        self.ref_depth = np.concatenate((np.arange(1,100,2), 
                                         np.arange(100,300,5), 
                                         np.arange(300, 1000, 50),
                                         np.arange(1000,4000,100)))
        
        # Subplot axes settings
        self.n_r = 2              # Number of subplot rows
        self.n_c = 9              # Number of subplot columns
        self.figsize = (6.5,5)      # Figure size
        self.sharey = True        # Align y axes
        self.sharex = False       # Align x axes
        self.subplot_padding = .4 # Amount of vertical and horizontal padding

        # Whole figure padding as % (left, bottom, right, top)
        self.fig_pad = (.09, .10, .15, .18)  
        self.max_depth = 150      # Maximum plot depth
        
        # Legend settings
        self.legend_str = ["CO9p2","CO7"]
        self.legend_fontsize = 8
        
        # Labels and Titles
        self.xlabelpos = (self.figsize[0]/2, 0)    # (x,y) position of xlabel
        
        self.label_fontsize = 6               # Fontsize of all labels
        self.label_fontweight = "normal"      # Fontweight to use for labels
                                              # and subtitles
        self.title_fontsize = 13              # Fontsize of title
        self.title_fontweight = "bold"        # Fontweight to use for title
        
        #%% SCRIPT: READ AND PLOT DATA
        
        # Read all datasets into list
        self.ds_list_DJF = [xr.open_dataset(dd) for dd in self.fn_list_DJF]
        self.ds_list_JJA = [xr.open_dataset(dd) for dd in self.fn_list_JJA]
        ds_list = self.ds_list_DJF
        self.n_ds = len(ds_list)
        self.n_reg = len(self.region_ind)
        
        print(f"Check region names specified are consistent with mask file")
        for i in range(self.n_reg):
            print(f"""Panel label:({self.region_names[i]}) matches data label:
                 ({ds_list[0].region_names.values[self.region_ind[i]]})""")
        
        # Loop over variable to plot
        for scalar in ["Temperature", "Salinity"]:
            for stat_type in ["MAE", "BIAS"]:
                self.plot_scalar_and_stat_type(scalar, stat_type)

    def plot_scalar_and_stat_type(self, scalar, stat_type):
        """
        Render plots for choice of scalar and statistical method.
        """

        if scalar == "Temperature": units = "($^{\circ}$C)"
        if scalar == "Salinity": units = "(-)"
    
        # Labels and Titles
        xlabel = "{0} (units)".format(stat_type)  # Xlabel string
    
        # Whole figure title
        fig_title = "Regional {0} || {1}".format(stat_type, scalar)
        
        if stat_type == "MAE":
            tmp_str = "mean_abs_diff"
        if stat_type == "STD":
            tmp_str = "std_diff"
        if stat_type == "BIAS":
            tmp_str = "mean_diff"
    
        var_name = "profile_{0}_{1}".format(tmp_str, scalar.lower())
    
        # Filename for the output
        fn_out = "FIGS/regional_{0}_{1}.pdf".format(var_name, self.run_name)
    
        # Create plot and flatten axis array
        f, axs = plt.subplots(self.n_r, self.n_c, figsize=self.figsize,
                              sharex=self.sharex, sharey=self.sharey)
    
        # Loop over regions
        for ii in range(self.n_reg):
            for row, season in enumerate(["DJF", "JJA"]):
                if ii >= self.n_reg:
                    axs[row,ii].axis('off')
                    continue
    
                # Get the index of this region
                index = self.region_ind[ii]
                
                # Loop over datasets and plot their variable
                p = []
                for pp in range(self.n_ds):
    
                    print(f"season:{season}")
                    if season == "DJF":
                      ds = self.ds_list_DJF[pp]
                    elif season == "JJA":
                      ds = self.ds_list_JJA[pp]
                    else:
                      print(f"Not expecting that season: {season}")
                    p.append(axs[row,ii].plot(ds[var_name][index][:100], 
                             self.ref_depth[:100])[0] )
    
                # Do some plot things
                axs[row,ii].set_title(f"{self.region_names[ii]}:\n{season}",
                                         fontsize=8)
                axs[row,ii].grid()
                axs[row,ii].set_ylim(0, self.max_depth)
    
                # set x lims
                if scalar == 'Salinity':
                    axs[row,ii].set_xlim(-0.1, 3.5)
                if scalar == 'Temperature':
                    axs[row,ii].set_xlim(-0.1, 4.0)
    
                # Plot fixed lines at 0 and mean depth
                if self.plot_zero_line:
                    axs[row,ii].plot([0,0], [0, self.max_depth], c='k',
                                        lw=0.5, ls='-')
                if self.plot_mean_depth:
                    axs[row,ii].plot([axs[row,ii].get_xlim()[0], 
                                      axs[row,ii].get_xlim()[1]], 
                                     [ds['bathymetry'][index],
                                      ds['bathymetry'][index]],
                                        color='k', ls='--')
    
                # Invert y axis
                axs[row,ii].invert_yaxis()
    
        # Make legend
        f.legend(p, self.legend_str, fontsize = self.legend_fontsize,
                 bbox_to_anchor=(1.0,0.5))
    
        # Set Figure title and label axes
        f.suptitle(fig_title, fontsize=self.title_fontsize,
                   fontweight=self.title_fontweight)
        f.text(0.5, 0.01, scalar + " " + units, ha='center')
        f.text(0.01, 0.5, 'Depth (m)', va='center', rotation=90)
    
        # Set x and y labels
        f.text(self.xlabelpos[0], self.xlabelpos[1], xlabel, 
             va = 'center', rotation='horizontal', 
             fontweight=self.label_fontweight, fontsize=self.label_fontsize)
    
        # Set tight figure layout
        f.tight_layout(w_pad=self.subplot_padding, h_pad=self.subplot_padding)
        f.subplots_adjust(left   = (self.fig_pad[0]), 
                          bottom = (self.fig_pad[1]),
                          right  = (1-self.fig_pad[2]),
                          top    = (1-self.fig_pad[3]))
    
        # Save plot maybe
        if self.save_plot: 
            f.savefig(fn_out)

    def plot_kattegat_norwegian_all_season(self):
        """
        Plot seasonal profiles of Kattegat and the Norwegian Trench.
        """

        # get data

        # initialise figure

        # plot

        # configure axes

        # add legend

        
if __name__ == "__main__":
    sp = seasonal_profiles()
    sp.plot_all_djf_jja()
