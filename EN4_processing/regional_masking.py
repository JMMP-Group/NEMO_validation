from PythonEnvCfg.config import config

import xarray as xr
import coast
from dask.diagnostics import ProgressBar
import numpy as np

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

        # set as global variable
        self.ds_option = ds_option

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
        
        mask_id = ['northern_north_sea','outer_shelf',
                   'eng_channel','nor_trench', 'kattegat', 'fsc',
                   'southern_north_sea', 'off_shelf', 'irish_sea' ]
        print("Size of names is ", len(mask_id))
        
        self.mask_xr = mm.make_mask_dataset(lon, lat, masks_list, mask_id)

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

        # tmp switch of dims for COAsT
        self.mask_xr = self.mask_xr.swap_dims({"region_names":"dim_mask"})

        # get mask indicies
        analysis = coast.ProfileAnalysis()
        mask_indices = analysis.determine_mask_indices(model_profile,
                                                            self.mask_xr)


        # make region names a coordinate dim
        self.mask_xr = self.mask_xr.swap_dims({"dim_mask":"region_names"})

        if self.ds_option == "bias": # difference between model and obs
            # get model profiles bathymetry - this seems clunky
            bathy_src = coast.Profile(config=self.cfg.fn_cfg_prof)
            fn_prof = self.cfg.dn_out + "profiles/{}_PRO_INDEX.nc".format(season)
            bathy_src.dataset = xr.open_dataset(fn_prof, chunks="auto")

            # Add bathymetry data to profiles before mask averaging
            model_profile.dataset['bathymetry'] = bathy_src.dataset.bathymetry

        # get stats per mask
        mask_stats = analysis.mask_stats(model_profile, mask_indices)

        # load stats for speed
        print ("loading mask stats")
        with ProgressBar():
            mask_stats = mask_stats.load()

        # make region names a coordinate dim
        mask_stats = mask_stats.swap_dims({"dim_mask":"region_names"})

        # extract regions
        model_profile_regions = []
        for region, region_mask in mask_indices.groupby("region_names"):

            mask_ind = region_mask.mask.astype(bool).squeeze()
            model_profile_region = model_profile.dataset.isel(id_dim=mask_ind)
            model_profile_region = model_profile_region.expand_dims(
                                      "region_names")
            model_profile_regions.append(model_profile_region)

        # combine extracted regions
        region_merged = xr.concat(model_profile_regions, dim='id_dim')

        return region_merged, mask_stats

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
        model_profile_seasons, mask_stats_seasons = [], []
        for season in seasons:
            print (f"Partitioning {season} by region")
            region_merged, stats = self.partition_profiles_by_region(
                                                 season=season)

            # add season coordinate
            region_merged = region_merged.assign_coords(
            {"season":("id_dim", np.full(region_merged.id_dim.shape, season))})

            # expand dims to include season
            stats = stats.expand_dims(season=[season])

            # add to list for conatenation
            model_profile_seasons.append(region_merged)
            mask_stats_seasons.append(stats)

        # combine by season
        model_profile_merged = xr.concat(model_profile_seasons, dim='id_dim')
        mask_stats_merged = xr.concat(mask_stats_seasons, dim='season')

        ## flatten
        #model_profile_flat = self.flatten_depth(model_profile_merged)

        # save
        with ProgressBar():
            print ("saving regional profiles")
            # region assigned profiles
            path = self.cfg.dn_out + self.fn_save + ".nc"
            model_profile_merged.to_netcdf(path)

            print ("saving regional stats")
            # regional stats
            path = self.cfg.dn_out + self.fn_save + "_stats.nc"
            mask_stats_merged.to_netcdf(path)

#    def create_masked_stats(self):
#        """
#        Use partitioned data to create masked statistics
#
#        None partition_by_region() should be run first
#        """
#
#     
#        # rename for COAsT
#        self.mask_indices = self.mask_indices.swap_dims(
#                                                    {"region_names":"dim_mask"})
#        self.model_profile_merged = self.model_profile_merged.swap_dims(
#                                                    {"region_names":"dim_mask"})
#
#        # Formatting for COAsT profile_analysis
#        # COAsT does not handle datetime for stats
#        # + error handleing does not catch > 1d arrays
#        time_type_list = ["datetime64[ns]","timedelta64[ns]"]
#        for var in self.model_profile_merged.data_vars:
#            if self.model_profile_merged[var].dtype in time_type_list:
#                self.model_profile_merged = self.model_profile_merged.drop_vars(
#                                                                       [var])
#
#        # masked stats requires COAsT objects
#        analysis = coast.ProfileAnalysis()
#        profile_object = coast.Profile(config=self.cfg.fn_cfg_prof)
#        profile_object.dataset = self.model_profile_merged
#
#
#        print (stats)
#
#        # save
#        with ProgressBar():


if __name__ == "__main__":
    ma = masking()
    ma.partition_by_region("bias")
#    ma.create_masked_stats()

    ma.partition_by_region("profiles")
#    ma.create_masked_stats()
    
