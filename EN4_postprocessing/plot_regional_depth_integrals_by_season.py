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
from EN4_processing.regional_masking import masking as r_mask
from EN4_postprocessing.plot_regional_mask import masking
import cartopy.crs as ccrs

matplotlib.rcParams.update({'font.size': 8})

class seasonal_depth_integral(object):
    '''
    Plotting collapsed measures of temperature and salinity biases per region.
    '''

    def __init__(self, case_num=1):

        self.case_num = case_num
        
        self.case_paths = [config.dn_out, config.comp_case["proc_data"]]
        self.models = [config.case, config.comp_case["case"]]

        # one or two model cases - trim
        self.case_paths = self.case_paths[:case_num]
        self.models = self.models[:case_num]

        self.ds_list = []
        for path in self.case_paths:
            fn = path + "/profiles/profile_bias_by_region_and_season_stats.nc"
            ds = xr.load_dataset(fn)

            # make region names indexable
            #ds = ds.swap_dims({"dim_mask":"region_names"})
            
            # append model list
            self.ds_list.append(ds)

        self.regions = ['northern_north_sea',
                   'outer_shelf',
                   'eng_channel',
                   'nor_trench',
                   'kattegat',
                   'southern_north_sea',
                   'irish_sea']
        self.clist = [plt.cm.tab10.colors[i] for i in [0,1,3,2,5,6,9]]
        self.region_names = ["N. North\nSea",
                        "Outer\nShelf",
                        "Eng.\nChannel",
                        "Nor.\nTrench", 
                        "Kattegat",
                        "S. North\nSea",
                        "Irish\nSea"]
    
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

        # alias var to globally accessible parameter
        self.var_str = scalar

        if scalar == "temperature": 
            x_label = "Temperature Bias ($^{\circ}$C)"
            y_max = 1.2
        if scalar == "salinity": 
            x_label = "Salinity Bias ($10^{-3}$)"
            y_max = 1.22
    
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

        # set bars
        self.render_bars(axs, self.da_list)

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

    def plot_regional_depth_integrals_bootstrapped_single_var(self,
                                                   scalar="temperature",
                                                   sample_size=1000):
        """
        Plot depth integrated differences between EN4 and NEMO.
        """

        # alias var to globally accessible parameter
        self.var_str = scalar

        if scalar == "temperature": 
            x_label = "Temperature Bias ($^{\circ}$C)"
            y_max = 1.4
        if scalar == "salinity": 
            x_label = "Salinity Bias ($10^{-3}$)"
            y_max = 2.5
    
        # initialise plot
        fig, axs = plt.subplots(1, figsize=(5.5,3.5))
        plt.subplots_adjust(top=0.98, right=0.98)
                              
        metric = "abs_diff"
        var = f"{metric}_{scalar}"
        fn = f"bootstrapped_{var}_bias_with_EN4_one_model_{sample_size}.nc"

        path_list = [config.dn_out+"profiles/" + fn,
                        config.comp_case["proc_data"] + "/profiles/"+ fn]
        path_list = path_list[:self.case_num]
        ds_list = [xr.load_dataset(dd) for dd in path_list]

        # select mean abs error for temperature or salinity  
        self.da_list = [] 
        for ds in ds_list:
            # select variable and depth average
            da = ds[f"{var}_quant_quant"]

            # list for all models
            self.da_list.append(da)

        # set bars
        self.render_bars(axs, self.da_list)

        axs.set_ylabel(x_label)
        axs.set_ylim(0,y_max) 
    
        # set transparent background
        fig.patch.set_alpha(0.0)

        # save
        model_str = ''
        for model in self.models:
            model_str += model + '_'
        save_name = "FIGS/" + model_str \
             + "bootstrapped_depth_integrated_regional_errors_by_season_cut_" \
                  + scalar + ".pdf"
        plt.savefig(save_name)

    def plot_regional_depth_integrals_bootstrapped_temp_salt(self):
        """
        Plot depth integrated differences between EN4 and NEMO for temperature
        and salinity.
        """

        # initialise plot
        #fig, axs = plt.subplots(2, figsize=(6.5,3.5))
        #plt.subplots_adjust(top=0.98, right=0.98)
        self.retrieve_record_stats()

        # initialise figure
        fig = plt.figure(figsize=(6.5,3.5))

        # initialise gridspec
        gs0 = gridspec.GridSpec(ncols=1, nrows=1)
        gs1 = gridspec.GridSpec(ncols=1, nrows=2)
    
        # set frame bounds
        gs0.update(top=0.90, bottom=0.1, left=0.04, right=0.35)
        gs1.update(top=0.98, bottom=0.1, left=0.44, right=0.99, 
                   hspace=0.1)

        # TODO check if mask exists
        mask_exists=False
        if mask_exists: # get masks
            mask_xr = xr.open_dataset(config.dn_out + "profiles/mask_xr.nc")
        else:
            mask_xr = r_mask().create_regional_mask()

        # initialise projection
        mid_lat = np.mean([44,mask_xr.latitude.max().values])
        mid_lon = np.mean([mask_xr.longitude.min().values,
                           mask_xr.longitude.max().values])
        mid_lon = -2.6
        proj_a=ccrs.EquidistantConic(central_latitude=mid_lat,
          standard_parallels=(44,mask_xr.latitude.max().values),
          central_longitude=mid_lon)

        # assign axes to lists
        axs = []
        axs.append(fig.add_subplot(gs0[0], projection=proj_a))
        for i in range(2):
            axs.append(fig.add_subplot(gs1[i]))

        # render mask
        self.render_regional_mask(axs[0], c_bar=True)

        # render bars
        legend_list = [True, False]
        for i, var in enumerate(["temperature","salinity"]):
            self.render_regional_depth_integrals_bootstrapped(axs[i+1], 
                                                  scalar=var,
                                                  add_legend=legend_list[i])

        # remove x-label from bar
        axs[1].set_xticklabels([])

        # add number of records to map
        proj=ccrs.PlateCarree()
        axins = axs[0].inset_axes([0.65, 0.02, 0.3, 0.2],
                         transform=axs[0].transAxes)
                          #xlim=inset_xlim, ylim=inset_ylim,
        n_records = self.retrieve_record_stats()
        x = np.arange(len(self.regions)) # the label locations
        for i, region in enumerate(self.regions):

            axins.bar(x[i], n_records[region], width=1.0, edgecolor=None,
                       linewidth=0,
                       facecolor=self.clist[i])
        axins.set_xticks([])
        axins.set_yticks([0,1e4])
        axins.set_yticklabels([0, 1], size=6)
        axins.spines[['right', 'top']].set_visible(False)
        axins.margins(x=0)
        axins.patch.set_alpha(0.0)
        axs[0].text(0.95, 0.22, "# of obs.\n" + r"(10$^4$)",
                    ha="right", va="top", size=6,
                    transform=axs[0].transAxes)



        # set transparent background
        fig.patch.set_alpha(0.0)

        # add panel labels
        axs[0].text(0.02, 0.98, "(a)", ha="left", va="top", size=6,
                    transform=axs[0].transAxes)
        axs[1].text(0.98, 0.96, "(b)", ha="right", va="top", size=6,
                    transform=axs[1].transAxes)
        axs[2].text(0.98, 0.96, "(c)", ha="right", va="top", size=6,
                    transform=axs[2].transAxes)

        # save
        model_str = ''
        for model in self.models:
            model_str += model + '_'
        save_name = "FIGS/" + model_str \
             + "bootstrapped_depth_integrated_regional_errors_by_season_cut_" \
                  + "temperature_and_salinity.pdf"
        plt.savefig(save_name)

    def retrieve_record_stats(self):
        """
        Get numbers of records for each region
        """

        path = config.dn_out + "profiles/profile_bias_by_region_and_season.nc"
        profiles = xr.open_dataset(path)
        n_records = {}
        for region in self.regions:
            region_prof = profiles.sel(region_names=region)
            nrec = region_prof.bathymetry.dropna("id_dim").sizes["id_dim"]
            n_records[region] = nrec

        return n_records

    def render_regional_mask(self, ax, c_bar=False):
        """
        Plot projected regional mask on existing ax.

        prm is masking class set in parent function
        """

        # initalise masking class
        prm = masking()

        # overwrite default regions in prm
        prm.regions = self.regions
        prm.clist = self.clist
        prm.region_names = self.region_names

        # get mask
        prm.get_mask()

        # get bathymetry - TODO: this information should be added to mask file
        prm.get_model_bathymetry()

        # projection setup
        proj=ccrs.PlateCarree()

        # initialise plot
        #fig, ax = plt.subplots(1, figsize=(5.5,3.5), subplot_kw=proj_dict)

        # render masks
        prm.render_regional_mask(ax, proj, c_bar=c_bar)

        # set axes labels
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")

    def render_regional_depth_integrals_bootstrapped(self, ax,
                                                     scalar="temperature",
                                                     sample_size=1000,
                                                     add_legend=False):
        """
        Render depth integrated differences between EN4 and NEMO on
        existing ax for a single variable.
        """

        self.var_str = scalar
        if scalar == "temperature": 
            x_label = "Temperature Bias ($^{\circ}$C)"
            y_max = 1.4
        if scalar == "salinity": 
            x_label = "Salinity Bias ($10^{-3}$)"
            y_max = 1.0
    
        metric = "abs_diff"
        var = f"{metric}_{scalar}"
        fn_in = f"bootstrapped_{var}_bias_with_EN4_one_model_{sample_size}.nc"

        path_list = [config.dn_out+"profiles/" + fn_in,
                        config.comp_case["proc_data"] + "/profiles/"+ fn_in]
        path_list = path_list[:self.case_num]
        ds_list = [xr.load_dataset(dd) for dd in path_list]

        # select mean abs error for temperature or salinity  
        self.da_list = [] 
        for ds in ds_list:
            # select variable and depth average
            da = ds[f"{var}_quant_quant"]

            # list for all models
            self.da_list.append(da)

        # set bars
        self.render_bars(ax, self.da_list, add_legend=add_legend)

        ax.set_ylabel(x_label)
        ax.set_ylim(0,y_max) 
    
    def get_bar_max_by_season(self, da_list):
        """ get index of maximum value between each model provided """

        # expand dims to have model id in order to merge
        da_list_n = []
        for i, da in enumerate(da_list):
            da_list_n.append(da.expand_dims(da_id=[i]))

        if len(da_list) > 1:
            # merge into dataset
            da = xr.merge(da_list_n)

            # get argmax
            da = da.argmax("da_id")

            # find which model has max for the first region per season
            max_da = da.sel(quantile=0.98,
                     region_names=self.regions[0]).to_dataarray().values[0]
        else:
            max_da = [0,0,0,0]
    
        return max_da


    def render_bars(self, ax, da_list, add_obs=False, add_legend=False):
        """ render season scatter coloured by region """


        for k, da in enumerate(da_list):
            # select regions
            da_list[k] = da.sel(region_names=self.regions)

        x = np.arange(len(self.regions)) # the label locations
        width = 0.2 / len(self.models)  # the width of the bars
         
        cmap = mcolors.ListedColormap(self.clist)
        seasons = ["DJF","MAM","JJA","SON"]

        # get index of max bars between models
        #if len(da_list) > 1:
        bar_max = self.get_bar_max_by_season(da_list)
        #else:
        #    bar_max = da.sel(quantile=0.98,
        #         region_names=self.regions[0]).values[0]

        # RDP - Too many loops, not readable...
        legend_data = []
        for k, da in enumerate(da_list):
            for j, season in enumerate(seasons):
                # select season
                bias = da.sel(season=season)
                for i, region in enumerate(self.regions):
                    offset = width * j * 1.2 * len(self.models) + (k * width)
                    bias_r = bias.sel(region_names=region)
                    if k > 0: 
                        alpha = 0.3
                    else:
                        alpha = 0.7 

                        if add_obs:
                            # add observational standard deviation
                            self.add_obs_std(ax, season, region,
                                             x[i] + offset, width * 2)

                    # render bars
                    bias_r_md = bias_r.sel(quantile=0.5).data
                    rect = ax.bar(x[i] + offset, bias_r_md, width,
                                  color=self.clist[i],
                                  alpha=alpha,
                                  align="edge")

                    # get data for legend
                    if (j==0) and (i==0):
                        legend_data.append(rect) 
                    
                    # render 96% confidence interval
                    lq = bias_r.sel(quantile=0.25).data
                    uq = bias_r.sel(quantile=0.75).data
                    vl = ax.vlines(x[i] + (width/2) + offset, lq, uq,
                             color=self.clist[i], transform=ax.transData,
                             lw=1)
                    vl.set(capstyle="round")

                    # add season labels
                    if i == 0 and k == bar_max[j]:
                        if k == 0:
                            # upper right
                            vertex = rect.patches[0].get_corners()[2]
                        else:
                            # upper left
                            vertex = rect.patches[0].get_corners()[3]
                        x_pos = vertex[0]
                        y_pos = vertex[1]

                        # set offset for season label
                        margin = ax.get_ylim()[1] * 0.05 * len(da_list)
                        
                        ax.text(x_pos, y_pos + margin, season, ha="center",
                                rotation=90, transform=ax.transData, size=6)

        # set x tick labels
        for i, name in enumerate(self.region_names):
            if name == "Kattegat":
                self.region_names[i] = "Katte."
        ax.set_xticks(x + (4*width*1.2*len(self.models) + width*k)/2,
                      self.region_names)

        if add_legend:
            ax.legend(legend_data, [config.case, config.comp_case["case"]],
                    loc='upper left', bbox_to_anchor=(0.02,0.96),
                            fontsize=6, borderaxespad=0, ncols=1)


    def add_obs_std(self, ax, season, region, x0, width):
        """
        add standard deviation of temperature/salinity to bar plot
        """

        # get data
        path = config.dn_out \
             + "masked_reductions/obs_season_merged_mask_std_mean.nc"
        obs_std = xr.open_dataset(path)[self.var_str + "_std"]

        # two standard deviations
        obs_std = obs_std * 2

        # set region_names as coordinate dimension
        #print (obs_std)
        #print (sdkjf)
        #obs_std = obs_std.set_index(dim_mask= "region_names")

        # select for region and season
        obs_std_region_season = obs_std.sel(seasons=season, region_names=region)

        # plot horizontal line
        ax.hlines(obs_std_region_season, x0, x0 + width,
                   transform=ax.transData, colors="k")

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
            mask_indices = self.get_mask_indicies(ds)
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
    
    def get_mask_indices(self, da, mask_exists=False):
        """
        Retrieve masked indicies

        Based on the assumption that mask_xr.nc has already been generated
        """

        # automate check for mask_xr.nc
        # TODO


        # define cfg files
        fn_cfg_nemo = config.fn_cfg_nemo
        fn_cfg_prof = config.fn_cfg_prof
        fn_dom_nemo = "%s%s"%(config.dn_dom, config.grid_nc)

        # get profile dataset
        obs_profiles = coast.Profile(config=fn_cfg_prof) 
        obs_profiles.dataset = da
        
        if mask_exists: # get masks
            mask_xr = xr.open_dataset(config.dn_out + "profiles/mask_xr.nc")
        else:
            mask_xr = r_mask().create_regional_mask()
        
        # get indices associated with each mask region
        analysis = coast.ProfileAnalysis()
        mask_indices = analysis.determine_mask_indices(obs_profiles,
                                                       mask_xr)
       
        return mask_indices.mask.astype(int)

    def format_to_box_plot(self, ds_var_quant, season, region):
        """ format data to conform with matplotlib bxp method """
        
        box = {
        'label' : f"{season} {region}",
        'whislo': ds_var_quant.sel(quantile=0.02).values, # 5th percentile
        'q1'    : ds_var_quant.sel(quantile=0.25).values, # 25th percentile
        'med'   : ds_var_quant.sel(quantile=0.50).values, # 50th percentile
        'q3'    : ds_var_quant.sel(quantile=0.75).values, # 75th percentile
        'whishi': ds_var_quant.sel(quantile=0.98).values, # 95th percentile
        "facecolor": "red"
        }

        return box

    def plot_angle_box_plot(self, scalar="temperature", sample_size=1000):
        """
        plot boxplot of bootstrapped statistics
        """

        metric = "abs_diff"
        var = f"{metric}_{scalar}"

        # get data
        fn = f"bootstrapped_{var}_bias_with_EN4_one_model_{sample_size}.nc"

        ds_path = config.dn_out + "profiles/" + fn
        ds_0 = xr.load_dataset(ds_path)

        ds_path = config.comp_case["proc_data"] + "profiles/" + fn
        ds_1 = xr.load_dataset(ds_path)

        da_0 = ds_0[f"{var}_quant_quant"]
        da_1 = ds_1[f"{var}_quant_quant"]

        # initiate plots
        fig, ax = plt.subplots(1, figsize=(6.5,4.5))
        plt.subplots_adjust(hspace=0.4,bottom=0.3,top=0.95)

        width=0.1
        boxes, pos = [], []
        #for i, (region_name, region) in enumerate(da_0.groupby("region_names")):
        #    for j, (season, subset) in enumerate(region.groupby("season")):
        
        for i, region in enumerate(da_0.region_names.data):
            for j, season in enumerate(da_0.season.data):
                for k, da in enumerate([da_0, da_1]):

                    subset = da.sel(region_names=region, season=season)
                    boxes.append(self.format_to_box_plot(subset.squeeze(),
                                                         season,
                                                         region))
                    pos.append( i +  (j * width * 2.2) + (k * width) )


        print (pos)
        # render
        bp = ax.bxp(boxes, positions=pos, widths=width,
                showfliers=False, patch_artist=True)
        
        clist = [plt.cm.tab10.colors[i] for i in range(9)]
        clist_rep = np.broadcast_to(clist, (8,9,3)).reshape(72,3, order="F")
        for i, (patch, color) in enumerate(zip(bp['boxes'],clist_rep)):
            patch.set_facecolor(color)
            if i % 2 == 1: # alpha based on model
                patch.set_alpha(0.4)

        # format y-axis
        ax.set_ylim([0,1.5])
        ax.set_ylabel(f"{scalar} bias")

        # format x-axis
        plt.xticks(rotation=60, ha='right')

        png_name = f"FIGS/{type}_bias_bootstrapped_box_plot.png"
        plt.savefig(png_name, dpi=600)

if __name__ == "__main__":
    sp = seasonal_depth_integral(case_num=2)
    sp.plot_regional_depth_integrals_bootstrapped_temp_salt()
