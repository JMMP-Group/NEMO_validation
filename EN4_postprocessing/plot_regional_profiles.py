from PythonEnvCfg.config import config
config = config() # initialise variables in python

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib
import string

matplotlib.rcParams.update({'font.size': 8})

class profiles(object):
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
        #co7_path = config.comp_case["proc_data"] + '/profiles/'
        self.fn_list = [config.dn_out+"profiles/" + fn,
                        config.comp_case["proc_data"] + "/profiles/"+ fn]

        #self.legend_str = [config.case,config.comp_case["case"]]
        #self.fn_list = [config.dn_out+"profiles/" + fn]
        self.legend_str = [config.case]
        self.n_ds = len(self.fn_list)
    
    def plot_all_regions(self, stat_type, stats_lim=True, save_plot=True):
        """
        plot profiles all regions and seasons
        
        2 rows and 7 columns with columns splitting region
         - row 1: temperature
         - row 2: salinity
        """

        # select regions
        self.region_id = ['northern_north_sea',
                          'outer_shelf',
                          'eng_channel',
                          'nor_trench',
                          'kattegat',
                          'southern_north_sea',
                          'irish_sea']

        # set region names
        self.region_names = ['N. North\nSea',
                             'Outer\nShelf',
                             'Eng.\nChannel',
                             'Nor.\nTrench',
                             'Kattegat',
                             'S. North\nSea',
                             'Irish\nSea']

        # initialise figure
        fig, axs = plt.subplots(2,7, figsize=(6.5,4.5))
        plt.subplots_adjust(left=0.1, right=0.98, top=0.93,
                            hspace=0.55, wspace=0.15, bottom=0.18)

        # get data
        self.ds_list_quant = [xr.load_dataset(dd.format("quants")) 
                        for dd in self.fn_list]
        self.ds_list_stats = [xr.load_dataset(dd.format("stats")) 
                        for dd in self.fn_list]

        # set line colours
        clist = [plt.cm.tab10.colors[i] for i in [0,1,3,2,5,6,9]]

        # choose metric
        if stat_type == "MAE":
            tmp_str = "mean_abs_diff"
            title_str = "Mean Abs. Err."
        if stat_type == "STD":
            tmp_str = "std_diff"
            title_str = "Std. Dev. Bias"
        if stat_type == "BIAS":
            tmp_str = "mean_diff"
            title_str = "Bias"

        # set units
        units = ["($^{\circ}$C)", "($10^3$)"]

        # Filename for the output
        fn_out = "FIGS/regional_profiles_{0}.pdf".format(config.case)

        ls=['-','--']
        lw=0.8
        leg_list = []
        for i, scalar in enumerate(["Temperature","Salinity"]):
            # choose variable
            var_name_quant = "{0}_{1}_quant_prof".format(tmp_str[5:],
                                                     scalar.lower())
            neg_lims, pos_lims = [], []
            for col, region in enumerate(self.region_id):
                axs[0,col].set_title(f"{self.region_names[col]}",
                                         fontsize=8, fontweight="bold")
                for mod in range(self.n_ds):
                    # plot MAE
                    ds_quant = self.ds_list_quant[mod].sel(season="ALL")
                    ds_stats = self.ds_list_stats[mod].sel(season="ALL")

                    # get region
                    var = ds_quant[var_name_quant].sel(region_names=region)

                    # restrict depth
                    region_bathy = ds_stats.profile_mean_bathymetry.sel(
                                   region_names=region)
                    var = var.where(var.depth < region_bathy)

                    lower_bound = var.sel(quantile=0.25)
                    mid_bound = var.sel(quantile=0.5)
                    upper_bound = var.sel(quantile=0.75)
                    
                    # get lims 
                    neg_lims.append(lower_bound.min())
                    pos_lims.append(upper_bound.max())

                    if mod==0:
                        axs[i,col].fill_betweenx(upper_bound.depth,
                                lower_bound, upper_bound,
                                        fc=clist[col], alpha=0.4)

                        p, = axs[i,col].plot(mid_bound, upper_bound.depth, lw=lw,
                                        c=clist[col], ls=ls[mod])
                        if (i == 0) and (col ==0): 
                            leg_list.append(p)
                    else:
                        axs[i,col].plot(upper_bound, upper_bound.depth, lw=lw,
                                        c='k', ls=ls[mod])
                        axs[i,col].plot(mid_bound, upper_bound.depth, lw=lw,
                                        c='k', ls='-')
                        p, = axs[i,col].plot(lower_bound, upper_bound.depth, 
                                        lw=lw, c='k', ls=ls[mod])
                        if (i == 0) and (col ==0): 
                            leg_list.append(p)

            # set value lims
            if stats_lim:
                glob_min = abs(np.quantile(np.array(neg_lims), 0.05))
                glob_max = abs(np.quantile(np.array(pos_lims), 0.95))
            else:
                glob_min = abs(min(neg_lims))
                glob_max = abs(max(pos_lims))
            bound = max(glob_min, glob_max) * 1.05
            for ax in axs[i].flatten():
                ax.set_xlim(-bound,bound)

            # set x-axis labels
            for ax in axs[i]:
                ax.set_xlabel(scalar + "\n" + title_str + "\n" + units[i])

        # general axes formatting
        for ax in axs.flatten():
            ax.axvline(0, lw=lw, c='grey')
            ax.set_ylim(0,100)
            ax.invert_yaxis()

        # blank out axes labels
        for ax in axs[:,1:].flatten():
            ax.set_yticklabels([])

        # set y-axis labels
        for ax in axs[:,0]:
            ax.set_ylabel("Depth (m)") 

        # add letters
        letters = ["({})".format(l) for l in 
                   list(string.ascii_lowercase)[:len(axs.flatten())]]
        print (letters)
        for i, ax in enumerate(axs.flatten()):
            ax.text(0.02, 0.02, letters[i], ha="left", va="bottom",
                    transform=ax.transAxes)

        # add legend
        axs[0,0].legend(leg_list, 
                 [config.case_label, config.comp_case["case_label"]],
                loc='lower right', bbox_to_anchor=(0.96,0.02),
                        fontsize=6, borderaxespad=0, ncols=1)

        # Save plot maybe
        if save_plot: 
            plt.savefig(fn_out)


if __name__ == "__main__":

    # set regions
    s_north_sea = {"region_id": "southern_north_sea",
                   "region_str": "S. North Sea"}
    irish_sea = {"region_id": "irish_sea",
                 "region_str": "Irish Sea"}

    # plot
    sp = profiles()
    #sp.plot_two_region_all_season(s_north_sea,
    #                              irish_sea,
    #                              scalar="temperature")
    #sp.plot_two_region_all_season(s_north_sea,
    #                              irish_sea,
    #                              scalar="salinity",
    #                              xlabel=r"$\overline{|\Delta S|}$ ($10^{-3}$)",
    #                              xmax=4.0)
    #sp.plot_all_djf_jja()
    sp.plot_all_regions(stat_type="BIAS")
