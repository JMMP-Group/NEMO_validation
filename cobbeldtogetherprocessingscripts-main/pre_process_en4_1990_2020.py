#!/usr/bin/env python3

import sys
# IF USING A DEVELOPMENT BRANCH OF COAST, ADD THE REPOSITORY TO PATH:
#sys.path.append('/data/users/fred/SET_UP_CONDA_ARTIFACTORY/COAST_SCIPY')
import coast
import xarray as xr
import numpy as np
import datetime
import pandas as pd

print('Modules loaded')

#fn_en4 = "/scratch/fred/EN4/EN.4.2.2.f.profiles.g10.*.nc"
fn_en4 = "/scratch/fred/EN4/EN.4.2.2.f.profiles.g10.*.nc"
fn_out = "/scratch/fred/EN4/SCIPY_processed_1990-2020.nc"
fn_cfg_prof = "/data/users/fred/coast_demo/config/example_en4_profiles.json"

profile = coast.Profile(fn_en4, multiple=True, config="./example_en4_profiles.json")

ind = profile.subset_indices_lonlat_box([-25.47, 16.25], [43, 64.5])[0]
profile = profile.isel(profile=ind)

new_profile = profile.process_en4()

new_profile.dataset.to_netcdf(fn_out)

