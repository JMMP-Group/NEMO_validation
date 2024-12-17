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

    def __init__(self):

        # paths
        self.fn_dom = cfg.dn_dom + cfg.grid_nc
        self.fn_path = cfg.dn_out + "surface_maps/"
        self.fn_comp_path = cfg.comp_case["proc_data"] + "surface_maps/"
        en4_nc = "/surface_maps/en4_gridded_surface_climatology.nc"
        m_nc = "/surface_maps/binned_model_surface_climatology.nc"
        self.en4_grid = cfg.dn_out + en4_nc
        self.model_binned = cfg.dn_out + m_nc

    def plot_bias_angle(self, var, depth_var=True):

        # TODO rename to "full depth" not "near full depth"
        fn=f"bootstrapped_{var}_bias_with_EN4.nc"

        # get data
        ds_path = cfg.dn_out + fn
        da = xr.open_dataset(ds_path)

        # initiate plots
        fig, axs = plt.subplots(9,4, figsize=(6.5,8.5))
        plt.subplots_adjust(hspace=0.4,bottom=0.05,top=0.95)
        colours = ["r","g","b","c"]
        for i, (region_name, region) in enumerate(da.groupby("region_names")):
            for j, (season, subset) in enumerate(region.groupby("season")):
                ax = axs[i,j]

                s = ax.plot(subset.bias_angle,
                            subset.temperature_frequency_quant.squeeze()[0],
                            c="slateblue", lw=0.8)
                s = ax.plot(subset.bias_angle,
                            subset.temperature_frequency_quant.squeeze()[1],
                            c="slateblue", ls="--", lw=0.8)
                s = ax.plot(subset.bias_angle,
                            subset.temperature_frequency_quant.squeeze()[2],
                            c="slateblue", lw=0.8)
                s = ax.plot(subset.bias_angle,
                           subset.temperature_frequency_mean.squeeze(),c='pink')
                ax.axvline(np.pi/4, c="orange")
                ax.axvline(-np.pi/4, c="orange")

                ax.set_title(f"{season} {region_name}")

        for ax in axs[:,0]:
            ax.set_ylabel("Freq.")

        for ax in axs[:,1:].flatten():
            ax.set_yticklabels([])

        for ax in axs.flatten():
            ax.set_xticks([-np.pi/2,-np.pi/4,0,np.pi/4,np.pi/2])

        for ax in axs[:-1,:].flatten():
            ax.set_xticklabels([])

        for ax in axs[-1,:]:
            ax.set_xticklabels(["$-\pi/2$",
                                "$-\pi/4$",
                                "0",
                                "$\pi/4$",
                                "$\pi/2$"])
            ax.set_xlabel(f"{var} bias")

        png_name = f"FIGS/full_depth_{var}_bias_angle_bootstrapped.png"
        plt.savefig(png_name, dpi=600)

bs =  bias_angle()
bs.plot_bias_angle("temperature")
