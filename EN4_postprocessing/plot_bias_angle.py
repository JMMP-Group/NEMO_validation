from PythonEnvCfg.config import config
cfg = config() # initialise variables in python

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xr
import numpy as np
import matplotlib

matplotlib.rcParams.update({'font.size': 8})

class bias_angle(object):
    """
    ploting routine for pointwise surface temperature and salinity maps
    """

    def __init__(self, var):

        # paths
        self.fn_dom = cfg.dn_dom + cfg.grid_nc
        self.fn_path = cfg.dn_out + "surface_maps/"
        self.fn_comp_path = cfg.comp_case["proc_data"] + "surface_maps/"
        en4_nc = "/surface_maps/en4_gridded_surface_climatology.nc"
        m_nc = "/surface_maps/binned_model_surface_climatology.nc"
        self.en4_grid = cfg.dn_out + en4_nc
        self.model_binned = cfg.dn_out + m_nc

        self.var = var
        self.metric = "abs_diff"
        self.type = f"{self.metric}_{self.var}"

    def plot_bias_angle_frequency(self, sample_size=1000):
        """
        Plot the bootstrapped bias angle frequency by season and region
        """

        # get data
        fn = f"bootstrapped_{self.type}_bias_with_EN4_{sample_size}.nc"
        ds_path = cfg.dn_out + fn
        da = xr.open_dataset(ds_path)

        # initiate plots
        fig, axs = plt.subplots(9,4, figsize=(6.5,8.5))
        plt.subplots_adjust(hspace=0.4,bottom=0.05,top=0.95)
        colours = ["r","g","b","c"]
        for i, (region_name, region) in enumerate(da.groupby("region_names")):
            for j, (season, subset) in enumerate(region.groupby("season")):
                ax = axs[i,j]
                var = f"{self.type}_frequency_quant"
                ax.fill_between(subset.bias_angle,
                subset[var].squeeze().sel(quantile=0.05),
                subset[var].squeeze().sel(quantile=0.95),
                                facecolor='slateblue', edgecolor=None,
                                alpha=0.5)
                ax.axvline(np.pi/4, c="orange", lw=0.8)
                med_var = \
            subset[f"{self.type}_angle_quant_quant"].squeeze().sel(quantile=0.5)
                mean_var = subset[f"{self.type}_angle_quant_mean"].squeeze()

                ax.axvline(med_var, c="red", lw=0.8)
                ax.axvline(mean_var, c="green", lw=0.8)

                ax.set_title(f"{season} {region_name}")

        for ax in axs[:,0]:
            ax.set_ylabel("Freq.")

        for ax in axs[:,1:].flatten():
            ax.set_yticklabels([])

        for ax in axs.flatten():
            ax.set_xlim([0,np.pi/2])
            ax.set_xticks([0,np.pi/4,np.pi/2])

        for ax in axs[:-1,:].flatten():
            ax.set_xticklabels([])

        for ax in axs[-1,:]:
            ax.set_xticklabels(["0",
                                "$\pi/4$",
                                "$\pi/2$"])
            ax.set_xlabel(f"{self.var} bias")

        png_name = f"FIGS/{self.type}_bias_angle_bootstrapped_histograms.png"
        plt.savefig(png_name, dpi=600)

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


    def plot_bias_angle_box_plot(self, sample_size=1000):
        """
        plot boxplot of bootstrapped statistics
        """

        # get data
        fn = f"bootstrapped_{self.type}_bias_with_EN4_{sample_size}.nc"
        ds_path = cfg.dn_out + fn
        ds = xr.open_dataset(ds_path, chunks=-1)
        da = ds.abs_diff_temperature_angle_quant_quant

        # initiate plots
        fig, ax = plt.subplots(1, figsize=(6.5,4.5))
        plt.subplots_adjust(hspace=0.4,bottom=0.3,top=0.95)

        boxes = []
        for i, (region_name, region) in enumerate(da.groupby("region_names")):
            for j, (season, subset) in enumerate(region.groupby("season")):

                boxes.append(self.format_to_box_plot(subset.squeeze(),
                                                     season,
                                                     region_name))

        ax.axhline(np.pi/4, c="slateblue", lw=0.8)

        # render
        bp = ax.bxp(boxes, showfliers=False,
                  patch_artist=True)

        
        clist = [plt.cm.tab10.colors[i] for i in range(9)]
        clist_rep = np.broadcast_to(clist, (4,9,3)).reshape(36,3, order="F")
        for patch,color in zip(bp['boxes'],clist_rep):
            patch.set_facecolor(color)

        # format y-axis
        ax.set_ylim([np.pi/16,5*np.pi/16])
        ax.set_yticks([np.pi/8,np.pi/4,3*np.pi/8])
        ax.set_yticklabels(["$\pi/8$","$\pi/4$","$3\pi/8$"])
        ax.set_ylabel(f"{self.var} bias")

        # format x-axis
        plt.xticks(rotation=60, ha='right')

        png_name = f"FIGS/{self.type}_bias_angle_bootstrapped_box_plot.png"
        plt.savefig(png_name, dpi=600)

if __name__ == "__main__":
    bs =  bias_angle("temperature")
    bs.plot_bias_angle_box_plot()
   # bs.plot_bias_angle_frequency()
