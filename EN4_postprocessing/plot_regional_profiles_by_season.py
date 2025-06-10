from PythonEnvCfg.config import config
config = config() # initialise variables in python

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib

matplotlib.rcParams.update({'font.size': 8})

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
        
        # Get two configurations
        fn = "profile_bias_by_region_and_season_{}.nc"
        co7_path = config.comp_case["proc_data"] + '/profiles/'
        self.fn_list = [config.dn_out+"profiles/" + fn,
                        config.comp_case["proc_data"] + "/profiles/"+ fn]

        self.legend_str = [config.case,config.comp_case["case"]]
        self.n_ds = len(self.fn_list)
    
    def plot_all_djf_jja(self):
        """
        Plot djf and jja for all locations.
        """

        # File settings
        self.run_name = "test"

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
        # TODO: THIS SHOULD BE ADDED DURING PROCESSSING...
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
        self.ds_list_stats = [xr.load_dataset(dd.format("stats")) 
                        for dd in self.fn_list]
        self.ds_list_quant = [xr.load_dataset(dd.format("quants")) 
                        for dd in self.fn_list]

        ds_list = self.ds_list_stats[0].sel(season="DJF")
        self.n_reg = len(self.region_ind)
        
        print(f"Check region names specified are consistent with mask file")
        for i in range(self.n_reg):
            print (i)
            print(f"""Panel label:({self.region_names[i]}) matches data label:
                 ({ds_list.region_names.values[self.region_ind[i]-1]})""")
        
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
    
        # TODO: variable naming needs to be consistent
        var_name_mean = "profile_{0}_{1}".format(tmp_str, scalar.lower())
        var_name_quant = "{0}_{1}_quant_prof".format(tmp_str[5:], scalar.lower())

    
        # Filename for the output
        fn_out = "FIGS/regional_{0}_{1}.pdf".format(var_name_mean, self.run_name)
    
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
                index = self.region_ind[ii] - 1
                
                # Loop over datasets and plot their variable
                p = []
                for pp in range(self.n_ds):
    
                    print(f"season:{season}")
                    if season in self.ds_list_stats[0].season:
                      ds_stats = self.ds_list_stats[pp].sel(season=season)
                      ds_quant = self.ds_list_quant[pp].sel(season=season)
                    else:
                      print(f"Not expecting that season: {season}")

                    # render mean profile
                    p.append(axs[row,ii].plot(
                             ds_stats[var_name_mean][index][:100], 
                             self.ref_depth[:100])[0] )


                    # get region
                    var = ds_quant[var_name_quant].isel(region_names=index)[:100]

                    # restrict depth
                    var = var.isel(z_dim=slice(None,100))

                    # render interquantile range
                    upper_bound = var.sel(quantile=0.95)
                    axs[row,ii].plot(upper_bound, upper_bound.depth, lw=0.8)
    
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
                                     [ds_stats['profile_mean_bathymetry'][index],
                                      ds_stats['profile_mean_bathymetry'][index]],
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

    def plot_two_region_all_season(
                         self,
                         r1_dict,
                         r2_dict,
                         scalar="temperature",
                         xmax=3.0, 
                         xlabel=r"$\overline{|\Delta \Theta|}$ ($^{\circ}$C)"):
        """
        Plot seasonal profiles of Kattegat and the Norwegian Trench.
        """
    
        # get data
        ds_list = [xr.load_dataset(dd) for dd in self.fn_list]

        ref_depth = np.concatenate((np.arange(1,100,2), 
                                    np.arange(100,300,5), 
                                    np.arange(300, 1000, 50),
                                    np.arange(1000,4000,100)))
        

        # initialise figure
        fig = plt.figure(figsize=(6.5,5.5))

        # initialise gridspec
        gs0 = gridspec.GridSpec(ncols=4, nrows=1)
        gs1 = gridspec.GridSpec(ncols=4, nrows=1)
    
        ## set frame bounds
        gs0.update(top=0.9, bottom=0.6, left=0.1, wspace=0.1, hspace=0.12,
                   right=0.85)
        gs1.update(top=0.4, bottom=0.1, left=0.1, wspace=0.1,right=0.85)

        # assign axes to lists
        axs = []
        for i in range(4):
            axs.append(fig.add_subplot(gs0[i]))
        for i in range(4):
            axs.append(fig.add_subplot(gs1[i]))

        ## plot
        def render(da, ax, season="DJF", region="", depth_lim=100):
            """
            render scalar on specified axis
            """

            # select season
            da = da.sel(season=season)

            # select region
            da = da.sel(region_names=region)

            # crude depth limit
            da_cut = da.isel(z_dim=slice(None,depth_lim))

            # render
            p, = ax.plot(da_cut["profile_mean_abs_diff_" + scalar].T,
                        ref_depth[:depth_lim])
            return p

        p_list = []
        for i, case in enumerate(self.legend_str):
            ds = ds_list[i]
            region_names = [r1_dict["region_id"],
                            r2_dict["region_id"]]
            for j, region in enumerate(region_names):
                render(ds, axs[0+4*(j-1)], "DJF", region)
                render(ds, axs[1+4*(j-1)], "MAM", region)
                render(ds, axs[2+4*(j-1)], "JJA", region)
                p0 = render(ds, axs[3+4*(j-1)], "SON", region)
                if j == 0:
                    p_list.append(p0)

        ## configure axes
        titles=["DJF","MAM","JJA","SON"]
        for i in range(4):
            axs[i].set_title(titles[i], size=8)
            axs[i+4].set_title(titles[i], size=8)
        for ax in axs[1:4] + axs[5:]:
            ax.set_yticklabels([])
        for ax in [axs[0], axs[4]]:
            ax.set_ylabel("Depth (m)") 
        for ax in axs:
            ax.set_ylim(0,200)
            ax.invert_yaxis()
            ax.set_xlim(0,xmax)
            ax.set_xlabel(xlabel)
        fig.legend(p_list, self.legend_str, fontsize=8,
                   bbox_to_anchor=(1.0,0.5))

        # add location annotation
        fig.text(0.475, 0.95, r1_dict["region_str"], 
             va = 'bottom', ha='center', rotation='horizontal', 
             fontsize=11)
        fig.text(0.475, 0.45, r2_dict["region_str"],
             va = 'bottom', ha='center', rotation='horizontal', 
             fontsize=11)

        # save
        save_path = "FIGS/" \
                    +"{}_{}_seasonal_mean_abs_diff_".format(
                    r1_dict["region_id"], r2_dict["region_id"]) \
                    + scalar + ".pdf"
        plt.savefig(save_path)
        
if __name__ == "__main__":

    # set regions
    s_north_sea = {"region_id": "southern_north_sea",
                   "region_str": "S. North Sea"}
    irish_sea = {"region_id": "irish_sea",
                 "region_str": "Irish Sea"}

    # plot
    sp = seasonal_profiles()
    #sp.plot_two_region_all_season(s_north_sea,
    #                              irish_sea,
    #                              scalar="temperature")
    #sp.plot_two_region_all_season(s_north_sea,
    #                              irish_sea,
    #                              scalar="salinity",
    #                              xlabel=r"$\overline{|\Delta S|}$ ($10^{-3}$)",
    #                              xmax=4.0)
    sp.plot_all_djf_jja()
