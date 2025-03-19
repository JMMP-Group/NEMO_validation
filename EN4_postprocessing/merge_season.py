from PythonEnvCfg.config import config
config = config() # initialise variables in python

import xarray as xr

def merge_mask_mean_seasons(path):
    """
    Merge seasonal regional errors in profiles of temperature and salinity

    TODO: this should be do at the point of file creation
    """

    fn_append = "_mask_means_daily.nc"
    seasons = []
    for season in ["DJF", "MAM", "JJA", "SON"]:
        # open dataset 
        ds_season = xr.open_dataset(
                             path + season + fn_append)

        # expand season dimension
        seasons.append(ds_season.expand_dims({"season": [season]}))

    # merge seasons
    merged = xr.merge(seasons)

    # save
    save_path = path + "season_merged_mask_means_daily.nc"
    merged.to_netcdf(save_path)
    print(f'file written to saved to {save_path}')
if __name__ == "__main__":
    co7_path = "/gws/nopw/j04/jmmp/CO9_AMM15_validation/co7/"
    co9_path = config.dn_out
    
    merge_mask_mean_seasons(co9_path)
    print("done")
