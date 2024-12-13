from PythonEnvCfg.config import config
cfg = config() # initialise variables in python

import xarray as xr
import numpy as np
from dask.diagnostics import ProgressBar

class bias_angle(object):
    """
    method for pointwise surface temperature and salinity histograms comparing
    bias across models referenced to EN4
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
        # TODO add season dimension
        fn="near_full_depth_EN4_bias_by_season_by_region.nc"
        
        # get data
        path = self.fn_path + fn
        self.da_0 = xr.open_dataset(path)["diff_" + self.var]

        path = self.fn_comp_path + fn
        self.da_1 = xr.open_dataset(path)["diff_" + self.var]

    def get_bias_angle_ds(self, depth_var=True):

        # flatten depth to make 2d (region,id_dim) array
        self.da_0 = self.flatten_across_depth(self.da_0)
        self.da_1 = self.flatten_across_depth(self.da_1)

        multiindex = ["season","time","latitude","longitude"]
        da_0 = self.da_0.set_index(id_dim=multiindex)#.drop_duplicates("id_dim")
        da_1 = self.da_1.set_index(id_dim=multiindex)#.drop_duplicates("id_dim")
        da_0, da_1 = xr.align(da_0, da_1)

        self.angle = np.arctan(da_0/da_1)
        self.angle.name = f"diff {self.var} angle"

    def flatten_across_depth(self, da):
        da = da.swap_dims({"z_dim":"depth"})
        da_depth = []
        for i, depth in enumerate(da.depth):
            da_depth.append(da.sel(depth=depth))
        da_1d = xr.concat(da_depth, dim="id_dim")

        return da_1d

    def get_bias_angle_hist(self, bootstrapped=True):
        """
        calculate 1-d histogram from bias angle data
        """

        hist_set = []
        for _, region in self.angle.groupby("region_names"):
            for season, subset in region.groupby("season"):
                if bootstrapped:
                    bootstrapped_hist = \
                                  self.get_bootstrapped_bias_angle_hist(subset,
                                                                        season)
                    hist_set.append(bootstrapped_hist)

        hist_set_ds = xr.merge(hist_set)

        fn = cfg.dn_out + f"bootstrapped_{self.var}_bias_with_EN4.nc"
        with ProgressBar():
            hist_set_ds.to_netcdf(fn)

    def get_bootstrapped_bias_angle_hist(self, ds, season, sample_size=1000):
        """
        Single instance bootstrap sample the bias angle and create histogram
        """

        bootstrap_set = []
        for i in range(sample_size):
            print (i)
            # randomised sampling
            obs_num = ds.sizes["id_dim"]
            random_sample = np.random.randint(obs_num, size=obs_num)
            ds_randomised = ds#ds.isel(id_dim=random_sample)

            # get histogram
            hist, bin_edges = np.histogram(ds_randomised, bins=100,
                                           range=[-np.pi/2,np.pi/2])
            bin_centres = (bin_edges[1:] + bin_edges[:-1]) / 2
            hist = np.expand_dims(hist, axis=[1,2])
  
            # assign to dataset
            bootstrapped_hist = xr.Dataset(
                          {f"{self.var}_frequency":(
                           ["bias_angle","region_names","season"],
                           hist)},
                          {"bias_angle":bin_centres,
                           "bias_angle_bin_edges": bin_edges,
                           "season":xr.DataArray([season], dims="season"),
                           "region_names":ds.region_names})
                           #["bias_angle","region_names","season","sample"],
                           #"sample":i})
            #bootstrapped_hist = bootstrapped_hist.expand_dims(dim={"sample":i})

            # append dataset to list
            bootstrap_set.append(bootstrapped_hist)

        # join boostrap hists as single dataset
        bootstrap_set_ds = xr.concat(bootstrap_set, dim="sample")

        # alias var name
        var_str = f"{self.var}_frequency"

        # get mean
        mean = bootstrap_set_ds.mean("sample")
        print (mean.to_dataarray())
        print (sdkf)
        mean = mean.rename({var_str:var_str + "_mean"})

        # get percentiles
        quant = bootstrap_set_ds.quantile([0.05, 0.5, 0.95], "sample")
        quant = quant.rename({var_str:var_str + "_quant"})
        
        # merge statistics
        bootstrap_statistics = xr.merge([mean, quant])

        return bootstrap_statistics

if __name__ == "__main__":
    ba = bias_angle("temperature")
    ba.get_bias_climatology()
    ba.get_bias_angle_ds()
    ba.get_bias_angle_hist()
