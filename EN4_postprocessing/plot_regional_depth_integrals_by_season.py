from PythonEnvCfg.config import config
config = config() # initialise variables in python

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib
import matplotlib.colors as mcolors
import coast
from coast import general_utils
from dask.diagnostics import ProgressBar

matplotlib.rcParams.update({'font.size': 8})

class seasonal_depth_integral(object):
    '''
    Plotting collapsed measures of temperature and salinity biases per region.
    '''

    def __init__(self, case_path, models):
        
        self.models = models

        self.ds_list = []
        for model in models:
            fn = model + "/profiles/season_merged_mask_means_daily.nc"
            ds = xr.open_dataset(case_path + fn)

            # make region names indexable
            ds = ds.swap_dims({"dim_mask":"region_names"})
            
            # append model list
            self.ds_list.append(ds)
    
    def depth_mean(self, da):

        # get depth at cell edges
        depth_w =  (da.depth[1:] + da.depth[:-1]) / 2
        
        # get cell thickness (e3t)
        e3t =  depth_w[1:] - depth_w[:-1]
        
        # extend to bottom cell
        e3t = xr.DataArray(np.concatenate(([e3t[0]], e3t, [e3t[-1]])),
                                           dims=("z_dim"))
            
        # depth mean
        return da.weighted(e3t).mean("z_dim")

        
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
        

        if scalar == "temperature": 
            x_label = "Temperature Bias ($^{\circ}$C)"
            y_max = 1.15
        if scalar == "salinity": 
            x_label = "Salinity Bias ($10^{-3}$)"
            y_max = 0.72
    
        # initialise plot
        fig, axs = plt.subplots(1, figsize=(5.5,3.5))
        plt.subplots_adjust(top=0.98, right=0.98)
                              
        # select mean abs error for temperature or salinity  
        self.da_list = []
        for ds in self.ds_list:
            # select variable and depth average
            da = self.depth_mean(ds["profile_mean_abs_diff_" + scalar])

            # list for all models
            self.da_list.append(da)

        # scatter
        self.render_bars(axs, self.da_list)

        # add obs profles std
        self.add_obs_std_to_bars(axs)

        axs.set_ylabel(x_label)
        axs.set_ylim(0,y_max) 
    
        # set transparent background
        fig.patch.set_alpha(0.0)

        # save
        model_str = ''
        for model in self.models:
            model_str += model + '_'
        save_name = "FIGS/" + model_str \
                  + "depth_integrated_regional_errors_by_season_cut_" \
                  + scalar + ".pdf"
        plt.savefig(save_name)

    def render_bars(self, ax, da_list):
        """ render season scatter coloured by region """

        regions = ['northern_north_sea',
                   'outer_shelf',
                   'eng_channel',
                   'nor_trench',
                   'kattegat',
                   'southern_north_sea',
                   'irish_sea']

        x = np.arange(len(regions)) # the label locations
        width = 0.2 / len(self.models)  # the width of the bars
         
        clist = [plt.cm.tab10.colors[i] for i in [0,1,3,2,5,6,9]]
        cmap = mcolors.ListedColormap(clist)
        seasons = ["DJF","MAM","JJA","SON"]
        # RDP - Too many loops, not readable...
        for k, da in enumerate(da_list):
            # select regions
            da = da.sel(region_names=regions)
            for j, season in enumerate(seasons):
                # select season
                bias = da.sel(season=season)
                for i, region in enumerate(regions):
                    #offset = width * ((k + 1) + j * 1.1)
                    offset = width * j * 1.2 * len(self.models) + (k * width)
                    bias_r = bias.sel(region_names=region)
                    if k > 0: 
                        alpha = 0.4
                    else:
                        alpha = 1
                    rect = ax.bar(x[i] + offset, bias_r, width, color=clist[i],
                                  alpha=alpha)
                    if i == 0 and k == 0:
                        ax.bar_label(rect, labels=[season], padding=3,
                                     rotation=90)
                        #x_pos = x[i] + offset + width
                        #right, top = rect.xy[0] + w, rect.xy[1] + h
                        #ax.text(top, x, labels=[season], padding=3,
                        #             rotation=90)

        region_names = ["N. North\nSea",
                        "Outer\nShelf",
                        "Eng.\nChannel",
                        "Nor.\nTrench", 
                        "Kattegat",
                        "S. North\nSea",
                        "Irish\nSea"]
        ax.set_xticks(x + (3*width*1.2*len(self.models) + width*k)/2,
                      region_names)

    def get_obs_std(self):
        """
        find the standard deviation of the interpolated obs profile and save
        TODO: this is at odds with regional_mean_by_season.py and these
        should be merged.
        The process has begun via addition of obs to extract_season.py 
        """

        def _preprocess(ds_month):
            """ drop broadcasting of depth variable """
            # TODO: this should be done in GEN_MOD_Dave_example_profile_vali...
            ds_month["depth"] = ds_month.depth.isel(id_dim=0)
            return ds_month

        # get observational profiles interpolated to uniform depths
        obs_path = config.dn_out + "profiles/interpolated_obs_*.nc"
        obs_profiles_all = xr.open_mfdataset(obs_path, combine='nested',
                                             concat_dim="id_dim",
                                             parallel=True,
                                             preprocess=_preprocess)

        # get standard deviation by season
        season_data = []
        for season, ds in obs_profiles_all.groupby("time.season"):
            # split by region
            mask_indices = self.get_mask_regions(ds)
            mask_data = ds.isel(id_dim=mask_indices)

            # standard devation for each depth
            mask_data_std = mask_data.std(["id_dim"], skipna=True)

            # depth-mean of standard devation
            mask_data_std_mean = self.depth_mean(mask_data_std)

            # set season dim
            mask_data_std_mean = mask_data_std_mean.expand_dims(
                                 dict(seasons=[season]))

            # append season list
            season_data.append(mask_data_std_mean)

        # join all season standard deviations
        std_seasons_regions = xr.concat(season_data, dim="seasons")

        # set variable names
        for var in std_seasons_regions.keys():
            std_seasons_regions = std_seasons_regions.rename({var:var+"_std"})

        # save
        with ProgressBar():
            std_seasons_regions.to_netcdf(config.dn_out
                     + "masked_reductions/obs_season_merged_mask_std_mean.nc")
    
    def get_mask_regions(self, da):
        """
        Retrieve masked regions

        Based on the assumption that mask_xr.nc has already been generated
        """

        # check for mask_xr.nc
        # TODO

        # define cfg files
        fn_cfg_nemo = config.fn_cfg_nemo
        fn_cfg_prof = config.fn_cfg_prof
        fn_dom_nemo = "%s%s"%(config.dn_dom, config.grid_nc)

        # get profile dataset
        obs_profiles = coast.Profile(config=fn_cfg_prof) 
        obs_profiles.dataset = da
        
        # get masks
        mask_xr = xr.open_dataset(config.dn_out + "profiles/mask_xr.nc")
        
        # get indices associated with each mask region
        analysis = coast.ProfileAnalysis()
        mask_indices = analysis.determine_mask_indices(obs_profiles,
                                                       mask_xr)
       
        return mask_indices.mask.astype(int)

if __name__ == "__main__":
    ds_path = "/gws/nopw/j04/jmmp/CO9_AMM15_validation/"
    sp = seasonal_depth_integral(ds_path, ["P2.0", "co7"])
    #sp.plot_regional_depth_integrals(scalar="temperature")
    sp.get_obs_std()
    #sp.plot_regional_depth_integrals(scalar="salinity")
