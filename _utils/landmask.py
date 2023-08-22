
def add_landmask(ds): 
    """
    Create a landmask array -- important for obs_operator. Calculated 
    from bottom_level.
 
    Input
    -----
    ds: xarray dataset with bottom level variable

    Output
    ------
    ds: xarray dataset with new landmask variable
    """

    # add landmask
    ds["landmask"] = ds.bottom_level == 0

    return ds
