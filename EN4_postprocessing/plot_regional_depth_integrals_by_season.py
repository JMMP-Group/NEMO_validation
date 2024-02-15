from PythonEnvCfg.config import config
config = config() # initialise variables in python

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib
import matplotlib.colors as mcolors

matplotlib.rcParams.update({'font.size': 8})

class seasonal_depth_integral(object):
    '''
    Plotting collapsed measures of temperature and salinity biases per region.
    '''

    def __init__(self, case_path):
        
        fn = "profiles/season_merged_mask_means_daily.nc"
        self.ds = xr.open_dataset(case_path + fn)

        # make region names indexable
        self.ds = self.ds.swap_dims({"dim_mask":"region_names"})
    
    def depth_mean(self):

        # Should match definition in EN4_processing: ref_depth
        # TODO: THIS SHOULD BE ADDED DURING PROCESSSING...
        ref_depth_t = np.concatenate((np.arange(1,100,2), 
                                    np.arange(100,300,5), 
                                    np.arange(300, 1000, 50),
                                    np.arange(1000,4000,100)))
 
        # interpolated e3t (crude estimate - w-pts = mid-pts between cells)
        e3t =  (ref_depth_t[1:] - ref_depth_t[:-1]) / 2
        
        # extend to bottom cell
        e3t = xr.DataArray(np.concatenate((e3t, [e3t[-1]])), dims=("z_dim"))
            
        # depth integral
        self.da = (self.da * e3t).sum("z_dim") / e3t.sum("z_dim")
        
    def plot_regional_depth_integrals(self, scalar="temperature"):
        """
        Plot depth integrated differences between EN4 and NEMO.
        """

        # Region indices (in analysis) to plot
        self.region_ind = [ 1, 7, 3, 2, 9, 5, 4, 6, 8]
        self.region_names = ['N. North\nSea','S. North\nSea',
                             'Eng.\nChannel','Outer\nShelf',
                             'Irish\nSea', 'Kattegat',
                             'Nor.\nTrench', 'FSC', 'Off-shelf']
        
        # select mean abs error for temperature or salinity  
        self.da = self.ds["profile_mean_abs_diff_" + scalar] 

        self.depth_mean()

        if scalar == "temperature": 
            x_label = "Temperature Bias ($^{\circ}$C)"
            y_max = 0.065
        if scalar == "salinity": 
            x_label = "Salinity Bias ($10^{-3}$)"
            y_max = 0.026
    
        # initialise plot
        fig, axs = plt.subplots(1, figsize=(5.5,3.5))
        plt.subplots_adjust(top=0.98, right=0.98)
                              
        # scatter
        self.render_bars(axs, self.da)

        axs.set_ylabel(x_label)
        axs.set_ylim(0,y_max) 
    
        # set transparent background
        fig.patch.set_alpha(0.0)

        # save
        save_name = "FIGS/depth_integrated_regional_errors_by_season_cut_{}.pdf".format(scalar)
        plt.savefig(save_name)

    def render_bars(self, ax, da):
        """ render season scatter coloured by region """

        regions = ['northern_north_sea',
                   'outer_shelf',
                   'eng_channel',
                   'nor_trench',
                   'kattegat',
                   'southern_north_sea',
                   'irish_sea']

        da = da.sel(region_names=regions)
        x = np.arange(len(regions)) # the label locations
        width = 0.2  # the width of the bars
         
        clist = [plt.cm.tab10.colors[i] for i in [0,1,3,2,5,6,9]]
        cmap = mcolors.ListedColormap(clist)
        for j, (season, bias) in enumerate(da.groupby("season")):
            for i, region in enumerate(regions):
                offset = (width * j * 1.1)
                bias_r = bias.sel(region_names=region)
                rect = ax.bar(x[i] + offset, bias_r, width, color=clist[i])
                if i == 0:
                    ax.bar_label(rect, labels=[season], padding=3, rotation=90)
        
        region_names = ['N. North\nSea',
                        'Outer\nShelf',
                        'Eng.\nChannel',
                        'Nor.\nTrench', 
                        'Kattegat',
                        'S. North\nSea',
                        'Irish\nSea']
        ax.set_xticks(x + (width*1.6), region_names)

if __name__ == "__main__":
    co7_path = '/gws/nopw/j04/jmmp/CO9_AMM15_validation/co7/'
    sp = seasonal_depth_integral(co7_path)
    sp.plot_regional_depth_integrals(scalar="temperature")
    sp.plot_regional_depth_integrals(scalar="salinity")
