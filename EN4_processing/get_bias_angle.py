from PythonEnvCfg.config import config
cfg = config() # initialise variables in python

import coast
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import numpy as np
import cmocean

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

    def get_bias_climatology(self):

        # TODO rename to "full depth" not "near full depth"
        fn="near_full_depth_EN4_bias_by_season_by_region.nc"
        
        # get data
        path = self.fn_path + fn
        self.da_0 = xr.open_dataset(path)["diff_" + self.var]

        path = self.fn_comp_path + fn
        self.da_1 = xr.open_dataset(path)["diff_" + self.var]

    def get_bias_angle_ds(self, depth_var=True):

        print (self.da_1)
        multiindex = ["season","time","latitude","longitude"]
        da_0 = self.da_0.set_index(id_dim=multiindex).drop_duplicates("id_dim")
        da_1 = self.da_1.set_index(id_dim=multiindex).drop_duplicates("id_dim")
        da_0, da_1 = xr.align(da_0, da_1)

        angle = np.arctan(da_0/da_1)
        angle.name = f"diff {var} angle"

    def get_bias_angle_hist(self):
        """
        calculate 1-d histogram from bias angle data
        """

        for j, region in enumerate(self.da_1.region_names.values):
            print ("region: ",  region)
            da_r0 = self.da_0.sel(region_names=region)
            da_r1 = self.da_1.sel(region_names=region)
            # get angle

            for i, season in enumerate(seasons):
                ax = axs[j,i]

                angle_season = angle.sel(season=season)

                if depth_var:
                    print (depth_var)
                    angle_season = angle_season.swap_dims({"z_dim":"depth"})

                    angle_season_depth = []
                    for i, depth in enumerate(angle.depth):
                        angle_season_depth.append(
                                                angle_season.sel(depth=depth))
                    angle_1d = xr.concat(angle_season_depth, dim="id_dim")


                    s = ax.hist(angle_1d, bins=100,
                                range=[-np.pi/2,np.pi/2])



    def get_bootstrapped_bias_angle_hist(self):
        """
        Single instance bootstrap sample the bias angle and create histogram
        """



if __name__ == "__main__":
    ba = bias_angle("temperature")
    ba.get_bias_climatology()
    ba.get_bias_angle_ds()
