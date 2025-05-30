from PythonEnvCfg.config import config

import xarray as xr
import coast
from dask.diagnostics import ProgressBar

class masking(object):
    """
    Mask plotting routines
    """

    def __init__(self):

        self.cfg = config() # initialise variables in python

        #%% File settings
        self.fn_cfg_nemo = self.cfg.fn_cfg_nemo
        self.fn_dom_nemo = self.cfg.dn_dom + self.cfg.grid_nc

        # open nemo lat/lon grid to define regions (as function of bathymetry)
        self.nemo = coast.Gridded(fn_domain=self.fn_dom_nemo,
                                  config=self.fn_cfg_nemo)

    def choose_data(self, ds_option):
        """ select profile dataset to process """

        if ds_option == "profiles": # model profiles
            self.fn_read = "profiles/{}_PRO_INDEX.nc"
            self.fn_save = "profiles/profiles_by_region_and_season"


        elif ds_option == "bias": # difference between model and obs
            self.fn_read = "profiles/{}_PRO_DIFF.nc"
            self.fn_save = "profiles/profile_bias_by_region_and_season"

        else:
            print ("Error: Option not implemented")

    def get_model_bathymetry(self):
        """ get nemo model bathymetry """
        
        # alias variables
        self.bath = self.nemo.dataset.bathymetry.values.squeeze()
        self.lon = self.nemo.dataset.longitude.values.squeeze()
        self.lat = self.nemo.dataset.latitude.values.squeeze()
        
        return self.lon, self.lat, self.bath

    def create_regional_mask(self):
        """
        Create regional mask

        This is replicated in regional_mean_by_season.py and
        merge_mean_crps.py and thus not needed here.
        TODO: make this a common function and remove from here?
        """
        
        # get model data
        lon, lat, bath = self.get_model_bathymetry()
       
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

        # make regions selectable
        self.mask_xr = self.mask_xr.swap_dims({"dim_mask":"region_names"})

    def partition_profiles_by_region(self, season=None):
        """ partition processed profiles by region """

        if season:
            fn_index = self.cfg.dn_out + self.fn_read.format(season)
            # get model profiles on EN4 grid
            model_profile = coast.Profile(config=self.cfg.fn_cfg_prof)
            model_profile.dataset = xr.open_dataset(fn_index, chunks="auto")
        else:
            # TODO add non seasonal handleing
            print ("non-seasonal arguments yet to be implemented")

        # create mask
        self.create_regional_mask()

        # get mask indicies
        analysis = coast.ProfileAnalysis()
        self.mask_indices = analysis.determine_mask_indices(model_profile,
                                                            self.mask_xr)

        # extract regions
        model_profile_regions = []
        for region, region_mask in self.mask_indices.groupby("region_names"):

            mask_ind = region_mask.mask.astype(bool).squeeze()
            model_profile_region = model_profile.dataset.isel(id_dim=mask_ind)
            model_profile_region = model_profile_region.expand_dims(
                                      "region_names")
            model_profile_regions.append(model_profile_region)

        # combine extracted regions
        region_merged = xr.concat(model_profile_regions, dim='id_dim')

        return region_merged

    def flatten_depth(self, da):

        # flatten depth to make 3d (season,region,id_dim) array
        da = da.swap_dims({"z_dim":"depth"})
        da_depth = []
        for i, depth in enumerate(da.depth):
            da_depth.append(da.sel(depth=depth))
        da_1d = xr.concat(da_depth, dim="id_dim")

        # make id_dim an multi-index dimension
        #multiindex = ["depth","time","latitude","longitude"]
        #da_1d = da_1d.set_index(id_dim=multiindex)

        # drop duplicate records (this removes union of regions!)
        #da_1d = da_1d.drop_duplicates("id_dim")

        #da_1d = da_1d.reset_index("id_dim")

        return da_1d

    def partition_by_region(self, ds="profiles"):
        """ 
        partition profile data by region and season over two nested loops
        """

        self.choose_data(ds)

        self.create_regional_mask()

        seasons = ["DJF","MAM","JJA","SON"]
        model_profile_seasons = []
        for season in seasons:
            print (f"Partitioning {season} by region")
            region_merged = self.partition_profiles_by_region(season=season)

            region_merged = region_merged.expand_dims(season=[season])

            model_profile_seasons.append(region_merged)

        # combine by season
        self.model_profile_merged = xr.concat(model_profile_seasons,
                                              dim='id_dim')

        ## flatten
        #model_profile_flat = self.flatten_depth(model_profile_merged)

        # save
        with ProgressBar():
            path = self.cfg.dn_out + self.fn_save + ".nc"
            self.model_profile_merged.to_netcdf(path)

    def create_masked_stats(self):
        """
        Use partitioned data to create masked statistics

        None partition_by_region() should be run first
        """

     
        # rename for COAsT
        self.mask_indices = self.mask_indices.swap_dims(
                                                    {"region_names":"dim_mask"})

        # Formatting for COAsT profile_analysis
        # COAsT does not handle datetime for stats
        # + error handleing does not catch > 1d arrays
        time_type_list = ["datetime64[ns]","timedelta64[ns]"]
        for var in self.model_profile_merged.data_vars:
            if self.model_profile_merged[var].dtype in time_type_list:
                self.model_profile_merged = self.model_profile_merged.drop_vars(
                                                                       [var])

        # masked stats requires COAsT objects
        analysis = coast.ProfileAnalysis()
        profile_object = coast.Profile(config=self.cfg.fn_cfg_prof)
        profile_object.dataset = self.model_profile_merged

        # get stats
        stats = analysis.mask_stats(profile_object, self.mask_indices)

        # save
        with ProgressBar():
            path = self.cfg.dn_out + self.fn_save + "_stats.nc"
            stats.to_netcdf(path)


if __name__ == "__main__":
    ma = masking()
    ma.partition_by_region("profiles")
    ma.create_masked_stats()
    
    ma.partition_by_region("bias")
    ma.create_masked_stats()
