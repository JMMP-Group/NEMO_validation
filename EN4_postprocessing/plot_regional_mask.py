from PythonEnvCfg.config import config

import xarray as xr
import matplotlib.pyplot as plt
import coast
import numpy as np
import cartopy.crs as ccrs
import matplotlib.colors as mcolors
import cartopy.feature as cfeature
import matplotlib
from dask.diagnostics import ProgressBar
from EN4_processing.regional_masking import masking as mask

matplotlib.rcParams.update({'font.size': 8})
plt.rcParams['figure.facecolor'] = 'black'

class masking(object):
    """
    Mask plotting routines
    """

    def __init__(self):

        self.cfg = config() # initialise variables in python

        #%% File settings
        self.fn_cfg_nemo = self.cfg.fn_cfg_nemo_t
        self.fn_dom_nemo = self.cfg.dn_dom + self.cfg.grid_nc

        # open nemo lat/lon grid to define regions (as function of bathymetry)
        self.nemo = coast.Gridded(fn_domain=self.fn_dom_nemo,
                                  config=self.fn_cfg_nemo)

        self.regions = ["northern_north_sea", "outer_shelf", "eng_channel",
                  "kattegat", "southern_north_sea", "off_shelf", "irish_sea"]
        self.region_names = ["N. North Sea", "Outer shelf","Eng. channel",
                        "Kattegat", "S. North Sea", "Off shelf", "Irish Sea" ]
        self.clist = [plt.cm.tab10.colors[i] for i in [0,1,3,5,6,8,9]]

    def quick_region_plot(self, mask: xr.Dataset):
        """
        Plot a map of masks in the MaskMaker object
        """
    
        n_mask = mask.dims["region_names"]
        offset = 10  # nonzero offset to make scaled-boolean-masks [0, >offset]
        for j in range(n_mask):
       	    tt = (j + offset) * mask["mask"].isel(region_names=j).squeeze()
            print (np.unique(tt))
       	    ff = tt.where(tt > 0).plot(x="longitude", y="latitude", 
                                       levels=range(offset,
                                       n_mask + offset + 1, 1), cmap='tab10',
                                       add_colorbar=False)

        cbar = plt.colorbar(ff)
        cbar.ax.get_yaxis().set_ticks([])
        region_names = ["N. North Sea", "Outer shelf","Eng. channel",
                        "Kattegat", "S. North Sea", "Off shelf", "Irish Sea" ]
        for j in range(0, n_mask, 1):
            cbar.ax.text(1 + 0.5,
                         offset + (j + 0.5),
                         region_names[j],
                         ha="left",
                         va="center",
                         color="k",
                         )
        cbar.ax.get_yaxis().labelpad = 15
        plt.title(None)

    def get_mask(self, mask_exists=False):
        """ open existing mask created by regional mean by season """

        if mask_exists: # get masks
            self.mask_xr = xr.open_dataset(config.dn_out + "profiles/mask_xr.nc")
        else:
            self.mask_xr = mask().create_regional_mask()
        
        # make regions selectable
        #self.mask_xr = self.mask_xr.swap_dims({"dim_mask":"region_names"})

    def get_model_bathymetry(self):
        """ get nemo model bathymetry """
        
        # alias variables
        self.bath = self.nemo.dataset.bathymetry.values.squeeze()
        self.lon = self.nemo.dataset.longitude.values.squeeze()
        self.lat = self.nemo.dataset.latitude.values.squeeze()
        
        return self.lon, self.lat, self.bath

    def render_regional_mask(self, ax, proj):
        """
        render projected mask to subplot panel
        """
    

        # subset regions
        ds = self.mask_xr.sel(region_names=self.regions)

        n_mask = ds.sizes["region_names"]
        offset = 10  # nonzero offset to make scaled-boolean-masks [0, >offset]
        cmap = mcolors.ListedColormap(self.clist)
        for j, region in enumerate(self.regions):
            if (region in ["nor_trench", "fsc"]) and \
              ("off_shelf" not in self.regions):
                continue
       	    tt = (j + 0.5) * ds.mask.sel(region_names=region).squeeze()
       	    mt = tt.where(tt > 0)
            ff = ax.contourf(ds.longitude, ds.latitude, mt,
                                       levels=range(0, n_mask+1), cmap=cmap,
                                       transform=proj)

        # add colour bar
        cbar = plt.colorbar(ff, pad=0.02, location="bottom",
                            orientation="horizontal")
        cbar.ax.get_xaxis().set_ticks([])
        for j in range(0, n_mask, 1):
            print (self.region_names[j].replace("\n"," "))
            cbar.ax.text(j+0.3, -0.3,
                         self.region_names[j].replace("\n"," "),
                         ha="left",
                         va="top",
                         rotation=-35,
                         color="k",
                         )

        if "fsc" in self.regions:
            c = plt.cm.tab10.colors[4]
            # add Faroe Shetland Channel
            ax.contour(self.mask_xr.longitude, self.mask_xr.latitude, 
                        self.mask_xr.mask.sel(region_names="fsc"),
                        colors=[c], transform=proj, linewidths=0.8)
            ax.annotate("FSC", (-5.5, 60.65), transform=proj, c=clist[1],
                         rotation=48, fontweight="bold")

        if ("nor_trench" in self.regions) and ("off_shelf" in self.regions):
            # add Norwegian Trench
            c = plt.cm.tab10.colors[2]
            ax.contour(self.mask_xr.longitude, self.mask_xr.latitude,
                        self.mask_xr.mask.sel(region_names="nor_trench"),
                        colors=[c], transform=proj, linewidths=0.8)
            ax.annotate("Nor. Trench", (3.7, 57.8), ha="center", va="center",
                          rotation=-26, transform=proj, c=c,
                          fontweight="bold", fontsize=6)

        # add land mask
        landmask = xr.where(self.nemo.dataset.bottom_level == 0, 1, np.nan)
        ax.contourf(landmask.longitude, landmask.latitude, landmask, 
                   colors=[plt.cm.Greys(0.2)],
                   transform=proj)
        
        # set extent - lon0 set according to presence of off shelf
        if "off_shelf" in self.regions: 
            lon0 = ds.longitude.isel(x_dim=0, y_dim=0).values
        else:
            lon0 = -12
        lon1 = 9.8
        ax.set_extent([lon0, lon1, 46, 62], ccrs.PlateCarree())

        # format gridlines
        lon_grid = [-20,-10, 0, 10]
        lat_grid = [45, 50, 55, 60, 65]
        gl = ax.gridlines(draw_labels=True, xlocs=lon_grid, ylocs=lat_grid,
                         color='k', alpha=0.5)
        gl.xpadding = 2
        gl.ypadding = 2
        gl.xlabel_style = {'size': 6}
        gl.ylabel_style = {'size': 6}
        gl.right_labels = False
        gl.bottom_labels = False
        plt.draw()
        
    def plot_regional_mask(self):
        """
        Plot projected regional mask
        """

        # get mask
        self.get_mask()

        # get bathymetry - TODO: this information should be added to mask file
        self.get_model_bathymetry()

        # projection setup
        proj=ccrs.PlateCarree()
        mid_lat = np.mean([44,self.mask_xr.latitude.max().values])
        mid_lon = np.mean([self.mask_xr.longitude.min().values,
                           self.mask_xr.longitude.max().values])
        mid_lon = -2.6
        proj_a=ccrs.EquidistantConic(central_latitude=mid_lat,
          standard_parallels=(44,self.mask_xr.latitude.max().values),
          central_longitude=mid_lon) 
        proj_dict = {"projection": proj_a}

        # initialise plot
        fig, ax = plt.subplots(1, figsize=(5.5,3.5), subplot_kw=proj_dict)
        plt.subplots_adjust(left=0.08, right=0.9, top=0.9,
                             bottom=0.1, wspace=0.05)
    
        # render masks
        self.render_regional_mask(ax, proj)

        # set axes labels
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")

        # set transparent background
        fig.patch.set_alpha(0.0)

        # save
        plt.savefig("FIGS/maskmaker_proj_plot", dpi=600)
    
    def plot_quick_regional_mask(self):
        """
        Plot unprojected regional mask.

        This does not use cartopy projections for speed.
        TODO: add projected version.
        """

        # get mask
        self.open_mask()

        # get bathymetry - TODO: this information should be added to mask file
        self.get_model_bathymetry()

        # plot
        #subset = ["N North Sea", "Outer shelf","Eng channel","Kattegat",
        #          "S North Sea", "Off shelf", "Irish Sea" ]
        subset = ["northern_north_sea", "outer_shelf", "eng_channel",
                  "kattegat", "southern_north_sea", "off_shelf", "irish_sea"]

        # render panel
        self.quick_region_plot(self.mask_xr.sel(region_names=subset))

        # add edge colour
        plt.contourf(self.lon,self.lat,self.bath, levels=(0,10), colors="w")
        plt.contour(self.lon,self.lat,self.bath, levels=(0,10), colors="k")

        # add Faroe Shetland Channel
        plt.contour(self.mask_xr.longitude, self.mask_xr.latitude, 
                    self.mask_xr.mask.sel(region_names="fsc"),
                    colors='r')
        plt.annotate("FSC", (-5.5, 60.75))

        # add Norwegian Trench
        plt.contour(self.mask_xr.longitude, self.mask_xr.latitude,
                    self.mask_xr.mask.sel(region_names="nor_trench"),
                    colors='r')
        plt.annotate("Nor. Trench", (4, 60.1), ha="center", va="center",
                      rotation=-60)

        # save
        plt.savefig("FIGS/maskmaker_quick_plot.png")

if __name__ == "__main__":
    mp = masking()
    #mp.plot_quick_regional_mask()
    mp.plot_regional_mask()
