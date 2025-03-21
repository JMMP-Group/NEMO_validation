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
        self.fn_path = cfg.dn_out + "profiles/"
        self.fn_comp_path = cfg.comp_case["proc_data"] + "profiles/"

        self.var = var
        self.metric = "abs_diff"
        self.type = f"{self.metric}_{self.var}"

    def get_bias_climatology(self):

        fn="profiles_by_region_and_season.nc"
        
        # get data
        path = self.fn_path + fn
        self.da_0 = xr.open_dataset(path, chunks="auto")[self.type]

        path = self.fn_comp_path + fn
        self.da_1 = xr.open_dataset(path, chunks="auto")[self.type]

    def get_bias_angle_ds(self, depth_var=True):

        
        region_group_0, region_group_1 = [], []
        for region in self.da_0.region_names:
            da_0_r = self.da_0.sel(region_names=region)
            da_1_r = self.da_1.sel(region_names=region)

            # make id_dim an multi-index dimension
            multiindex = ["time","latitude","longitude"]
            da_0_r = da_0_r.set_index(id_dim=multiindex)
            da_1_r = da_1_r.set_index(id_dim=multiindex)

            # drop duplicate records - why are there duplicates?
            da_0_r = da_0_r.dropna("id_dim", how="all")
            da_1_r = da_1_r.dropna("id_dim", how="all")

            # align
            da_0_r, da_1_r = xr.align(da_0_r, da_1_r)

            da_0_r = da_0_r.expand_dims("region_names")
            da_1_r = da_1_r.expand_dims("region_names")

            region_group_0.append(da_0_r)
            region_group_1.append(da_1_r)

        # regroup regions
        da_0 = xr.concat(region_group_0, dim="id_dim")
        da_1 = xr.concat(region_group_1, dim="id_dim")

        # get angle
        self.angle = np.arctan(da_0/da_1)
        self.angle.name = f"{self.type}_angle"

        self.angle = self.angle.reset_index("id_dim")

        # save
        fn = cfg.dn_out + f"{self.type}_bias_with_EN4.nc"
        with ProgressBar():
            self.angle.to_netcdf(fn)

    def get_bias_angle_hist(self, bootstrapped=True):
        """
        calculate 1-d histogram from bias angle data
        """

        sample_size=1000

        # get data
        fn = cfg.dn_out + f"{self.type}_bias_with_EN4.nc"
        angle = xr.open_dataarray(fn, chunks="auto")

        hist_set = []
        for j, (_, region) in enumerate(angle.groupby("region_names")):
            for i, (season, subset) in enumerate(region.groupby("season")):
                print (i, j)
                if bootstrapped:
                    bootstrapped_hist = \
                                  self.get_bootstrapped_bias_angle_hist(subset,
                                                                    sample_size)
                    hist_set.append(bootstrapped_hist)

        hist_set_ds = xr.merge(hist_set)

        # set attrs
        hist_set_ds = hist_set_ds.assign_attrs(
                      {"bootstrap_sample_size":sample_size})

        fn = cfg.dn_out + f"bootstrapped_{self.type}_bias_with_EN4_{sample_size}.nc"
        with ProgressBar():
            hist_set_ds.to_netcdf(fn)

    def get_bootstrapped_bias_angle_hist(self, ds, sample_size=1000):
        """
        Single instance bootstrap sample the bias angle and create histogram
        """

        ds = ds.load()

        bootstrap_set = []
        for i in range(sample_size):

            # remove nans introduced during season-region broadcasting
            ds = ds.dropna("id_dim", how="all")

            # randomised sampling
            obs_num = ds.sizes["id_dim"]
            random_sample = np.random.randint(obs_num, size=obs_num)
            ds_randomised = ds.isel(id_dim=random_sample)

            # get histogram
            hist, bin_edges = np.histogram(ds_randomised, bins=100,
                                           range=[-np.pi/2,np.pi/2])
            bin_centres = (bin_edges[1:] + bin_edges[:-1]) / 2
            hist = np.expand_dims(hist, axis=[1,2])

            quant = ds_randomised.quantile(0.5, ["id_dim","z_dim"])
            quant.name = f"{self.type}_angle_quant"
            quant = quant.rename({"quantile":"sample_quant"})
            mean = ds_randomised.mean(["id_dim","z_dim"])
            mean.name = f"{self.type}_angle_mean"

            # merge seasons
            agg_stats = xr.merge([mean,quant])
  
            # assign to dataset
            bootstrapped_hist = xr.Dataset(
                          {
                           f"{self.type}_frequency":(
                           ["bias_angle","region_names","season"],
                           hist)
                           },
                          {"bias_angle":bin_centres,
                           "bias_angle_bin_edges": bin_edges,
                           "season":ds.season,
                           "region_names":ds.region_names})

            # merge stats with hist
            bootstrapped_vars = xr.merge([bootstrapped_hist, agg_stats])

            # append dataset to list
            bootstrap_set.append(bootstrapped_vars)

        # join boostrap hists as single dataset
        bootstrap_set_ds = xr.concat(bootstrap_set, dim="sample").load()

        # alias var name
        var_str = f"{self.type}_frequency"

        # get mean
        mean = bootstrap_set_ds.mean("sample")
        mean = mean.rename({var_str:var_str + "_mean",
           f"{self.type}_angle_mean":f"{self.type}_angle_mean_mean",
           f"{self.type}_angle_quant":f"{self.type}_angle_quant_mean"})

        # get percentiles
        quantiles = [0.02,0.05,0.25,0.5,0.75,0.95,0.98]
        quant = bootstrap_set_ds.quantile(quantiles,"sample")
        quant = quant.rename({var_str:var_str + "_quant",
           f"{self.type}_angle_mean":f"{self.type}_angle_mean_quant",
           f"{self.type}_angle_quant":f"{self.type}_angle_quant_quant"})
        
        # merge statistics
        bootstrap_statistics = xr.merge([mean, quant])

        return bootstrap_statistics

if __name__ == "__main__":
    ba = bias_angle("temperature")
    ba.get_bias_climatology()
    #ba.get_bias_angle_ds()
    ba.get_bias_angle_hist()
