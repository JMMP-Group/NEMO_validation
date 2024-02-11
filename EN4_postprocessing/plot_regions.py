from PythonEnvCfg.config import config
config = config() # initialise variables in python

import xarray as xr
import matplotlib.pyplot as plt
import coast

class mask_plotting(object):
    """
    Mask plotting routines
    """

    def __init__(self):
        #%% File settings
        fn_cfg_nemo = config.fn_cfg_nemo
        fn_dom_nemo = config.dn_dom + config.grid_nc

    def region_plot(self, mask: xr.Dataset):
        """
        Plot a map of masks in the MaskMaker object
        Add labels
        """
    
        n_mask = mask.dims["dim_mask"]
        offset = 10  # nonzero offset to make scaled-boolean-masks [0, >offset]
        for j in range(0, n_mask, 1):
       	    tt = (j + offset) * mask["mask"].isel(dim_mask=j).squeeze()
       	    ff = tt.where(tt > 0).plot(x="longitude", y="latitude", 
                                       levels=range(offset,
                                       n_mask + offset + 1, 1), cmap='tab10',
                                       add_colorbar=False)

        cbar = plt.colorbar(ff)
        cbar.ax.get_yaxis().set_ticks([])
        for j in range(0, n_mask, 1):
            cbar.ax.text(1 + 0.5,
                         offset + (j + 0.5),
                         mask["region_names"].isel(dim_mask=j).values,
                         ha="left",
                         va="center",
                         color="k",
                         )
        cbar.ax.get_yaxis().labelpad = 15
        plt.title(None)

    def open_mask(self):
        """ open existing mask created by regional mean by season """

        path = config.dn_out + "mask_xr.nc"
        self.regional_mask = xr.open_dataset(path, chunks=-1)

    def get_model_bathymetry(self):
        """ get nemo model bathymetry """
        
        # load nemo lat/lon grid to define regions (as function of bathymetry)
        nemo = coast.Gridded(fn_domain=fn_dom_nemo, config=fn_cfg_nemo)
        
        # alias variables
        self.bath = nemo.dataset.bathymetry.values.squeeze()
        self.lon = nemo.dataset.longitude.values.squeeze()
        self.lat = nemo.dataset.latitude.values.squeeze()

    def create_regional_mask(self, overlap=False):
        """
        Create regional mask

        This is replicated in regional_mean_by_season.py and
        merge_mean_crps.py and thus not needed here.
        TODO: make this a common function and remove from here?
        """
        
        # get model data
        self.get_model_bathymetry(self)
       
        # Define Regional Masks
        mm = coast.MaskMaker()
        masks_list = []

        # Add regional masks
        masks_list.append(mm.region_def_nws_north_north_sea(lon, lat, bath))
        masks_list.append(mm.region_def_nws_outer_shelf(lon, lat, bath))
        masks_list.append(mm.region_def_nws_english_channel(lon, lat, bath))
        masks_list.append(mm.region_def_nws_norwegian_trench(lon, lat, bath))
        masks_list.append(mm.region_def_nws_kattegat(lon, lat, bath))
        masks_list.append(mm.region_def_nws_fsc(lon, lat, bath))
        masks_list.append(mm.region_def_nws_south_north_sea(lon, lat, bath))
        masks_list.append(mm.region_def_nws_off_shelf(lon, lat, bath))
        masks_list.append(mm.region_def_nws_irish_sea(lon, lat, bath))
        
        masks_names = ["N North Sea", "Outer shelf","Eng channel",
                       "Norwegian Trench", "Kattegat", "FSC",
                       "S North Sea", "Off shelf", "Irish Sea" ]
        print("Size of names is ", len(masks_names[:]))
        
        self.mask_xr = mm.make_mask_dataset(lon, lat, masks_list, masks_names)

    def plot_regional_mask(self):
        """ plot regional mask """

        # get mask
        self.open_mask()

        # get bathymetry - TODO: this information should be added to mask file
        self.get_model_bathymetry()

        # plot
        region_plot(mask_xr.drop(["Norwegian Trench", "FSC"]))
        plt.contourf(lon,lat,bath, levels=(0,10), colors="w")
        plt.contour(lon,lat,bath, levels=(0,10), colors="k")

        # add Faroe Shetland Channel
        plt.contour(mask_xr.longitude, mask_xr.latitude, 
                    mask_xr.mask.isel(dim_mask=7),colors='r')
        plt.annotate("FSC", (-5.5, 60.75))

        # add Norwegian Trench
        plt.contour(mask_xr.longitude, mask_xr.latitude,
                    mask_xr.mask.isel(dim_mask=8),colors='r')
        plt.annotate("Nor. Trench", (4, 60.1), ha="center", va="center",
                      rotation=-60)

        # save
        plt.savefig("FIGS/maskmaker_quick_plot.png")

if __name__ == "__main__":
    mp = mask_plotting()
    mp.plot_regional_mask()
