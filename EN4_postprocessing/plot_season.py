from EN4_postprocessing.plot_regional_profiles_by_season import seasonal_profiles
from EN4_postprocessing.plot_regional_mask import masking
from EN4_postprocessing.plot_regional_depth_integrals_by_season import seasonal_depth_integrals

"""
Plot seasonal profiles consistent with Byrne et al. (2023) + extra
"""

# plot profiles for summer and winter
sp = seasonal_profiles()
sp.plot_all_djf_jja()

# plot regional mask
mp = masking()
mp.plot_regional_mask()

# depth integrated bias plotting
sp = seasonal_depth_integral()
sp.get_obs_std()
sp.plot_regional_depth_integrals(scalar="temperature")
sp.plot_regional_depth_integrals(scalar="salinity")
