'''
For plotting the analysis regions.
Uses COAsT example data
'''



import xarray as xr
import matplotlib.pyplot as plt
import sys

coast_repo = "/Users/jelt/GitHub/COAsT"
sys.path.append(coast_repo)
import coast


def region_plot(mask: xr.Dataset):
	"""
    Plot a map of masks in the MaskMaker object
    Add labels
    """

	n_mask = mask.dims["dim_mask"]
	offset = 10  # nonzero offset to make scaled-boolean-masks [0, >offset]
	for j in range(0, n_mask, 1):
		tt = (j + offset) * mask["mask"].isel(dim_mask=j).squeeze()
		ff = tt.where(tt > 0).plot(
			x="longitude", y="latitude", levels=range(offset, n_mask + offset + 1, 1), cmap='tab10', add_colorbar=False
		)

	cbar = plt.colorbar(ff)
	cbar.ax.get_yaxis().set_ticks([])
	for j in range(0, n_mask, 1):
		cbar.ax.text(
			1 + 0.5,
			offset + (j + 0.5),
			mask["region_names"].isel(dim_mask=j).values,
			ha="left",
			va="center",
			color="k",
		)
	cbar.ax.get_yaxis().labelpad = 15
	plt.title(None)

#%% File settings
fn_dom_nemo = "/Users/jelt/GitHub/COAsT/example_files/coast_example_nemo_domain.nc"
fn_cfg_nemo = "/Users/jelt/GitHub/COAsT/config/example_nemo_grid_t.json"

print('Doing regional analysis..')

# load nemo lat/lon grid to define regions (as function of bathymetry)
nemo = coast.Gridded(fn_domain=fn_dom_nemo, config=fn_cfg_nemo)

bath = nemo.dataset.bathymetry.values.squeeze()
lon = nemo.dataset.longitude.values.squeeze()
lat = nemo.dataset.latitude.values.squeeze()


# Define Regional Masks

mm = coast.MaskMaker()
masks_list = []
# Add regional mask for whole domain
#masks_list.append(np.ones(lon.shape))  # 0
masks_list.append(mm.region_def_nws_north_north_sea(lon, lat, bath))  # 1
masks_list.append(mm.region_def_nws_outer_shelf(lon, lat, bath))  # 2
masks_list.append(mm.region_def_nws_english_channel(lon, lat, bath))  # 3
#masks_list.append(mm.region_def_nws_norwegian_trench(lon, lat, bath))  # 4
masks_list.append(mm.region_def_nws_kattegat(lon, lat, bath))  # 5
#masks_list.append(mm.region_def_nws_fsc(lon, lat, bath))  # 6
masks_list.append(mm.region_def_nws_south_north_sea(lon, lat, bath))  # 7
masks_list.append(mm.region_def_nws_off_shelf(lon, lat, bath))  # 8
masks_list.append(mm.region_def_nws_irish_sea(lon, lat, bath))  # 9

#masks_names = ['whole_domain', 'northern_north_sea', 'outer_shelf', 'eng_channel', 'nor_trench',
#			   'kattegat', 'fsc', 'southern_north_sea', 'off_shelf', 'irish_sea']
masks_names = ['N North Sea','Outer shelf','Eng channel','Kattegat',
	     'S North Sea', 'Off shelf', 'Irish Sea' ]
print("Size of names is ", len(masks_names[:]))

mask_xr = mm.make_mask_dataset(lon, lat, masks_list, masks_names)


region_plot(mask_xr)
plt.contourf(lon,lat,bath, levels=(0,10), colors="w")
plt.contour(lon,lat,bath, levels=(0,10), colors="k")

## Append regions for contouring. These overlap shaded regions
masks_list.append(mm.region_def_nws_norwegian_trench(lon, lat, bath))  # 5
masks_list.append(mm.region_def_nws_fsc(lon, lat, bath))  # 6
masks_names = [ 'N North Sea', 'Outer shelf', 'Eng Channel', 'Kattegat',
			    'S North_Sea', 'Off shelf', 'Irish Sea', 'Nor Trench', 'FSC']
mask_xr = mm.make_mask_dataset(lon, lat, masks_list, masks_names)

plt.contour(mask_xr.longitude, mask_xr.latitude, mask_xr.mask.isel(dim_mask=7),colors='r')
plt.annotate("FSC", (-5.5, 60.75))
plt.contour(mask_xr.longitude, mask_xr.latitude, mask_xr.mask.isel(dim_mask=8),colors='r')
plt.annotate("Nor. Trench", (4, 60.1), ha="center", va="center", rotation=-60)

plt.savefig("FIGS/maskmaker_quick_plot.png")
plt.close("all")
