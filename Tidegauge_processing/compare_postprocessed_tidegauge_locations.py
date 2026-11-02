# Match CO9 port analysis to match CO7 analysis locations.
"""
There are 51 ports in CO7 and 56 ports in the CO9 (P1.5c) analysis
There are 61 ports in the obs dataset.

The analysis on CO7 was done some years ago and the parent data is recoverable but very slow to do so. Can we use the existing CO7 analysis?

After stripping out duplicates there are still 51 ports in the CO7 analysis, and 53 in the CO9 analysis.
The extra ports in the reduced CO9 can be identified and removed.

Method:
. demonstrate the co7 ports have no duplicates
. for each co7 port find a match in co9
. save the reduced set of co9 ports.

jelt 18 Dec 2024
"""

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

dir = "/gws/nopw/j04/jmmp/CO9_AMM15_validation/"
fn_co7=dir + "co7/tg_analysis/ssh_hourly_analyse_co7.nc"
fn_co9=dir + "P1.5c/tg_analysis/ssh_hourly_analyse_P1.5c.nc"
fn_obs="/gws/nopw/j04/jmmp/CO9_AMM15/obs/tg_amm15.nc"

dir = "/Users/jelt/Downloads/"
fn_co7 = dir + "ssh_hourly_analyse_co7.nc"
fn_co9 = dir + "ssh_hourly_analyse_P1.5c.nc"
#fn_co9 = dir + "ssh_hourly_analyse_P1.5c_matchCO7.nc"  # Use this file to check it worked
fn_obs = dir + "tg_amm15.nc"

########################
# load data

ds7 = xr.open_dataset(fn_co7)
ds9 = xr.open_dataset(fn_co9)
dsobs = xr.open_dataset(fn_obs)  # actually has a different data structure. Is only used here to check number of ports
#dsobs

print(f" co7 locations: {len(ds7.id_dim)}\n co9 locations: {len(ds9.id_dim)}\n obs locations: {len(dsobs.port)}")



########################
# Are the locations in CO7 unique?
#  Use complex numbers to pair coords. Find unique pairs
def unique_coords(list_cor):
    # list_cor   [[4190091.4195999987, 7410226.618699998], 
    #[4190033.2124999985, 7410220.0823] ...]
    coords=[ x+1j*y for (x,y) in list_cor] # using complex; a way for grouping
    uniques, ind, counts=np.unique(coords, return_index=True, return_counts=True)
    res=[ [x.real,x.imag] for x in uniques ] # ungroup
    xx = [ x.real for x in uniques ] 
    yy = [ x.imag for x in uniques ] 
    return res , xx ,yy, ind

#tt = list(ds9.id_dim)
#tt.remove(52)
#tt.remove(30)

list_coords_ds9 = [ [float(ds9.longitude[i].values), float(ds9.latitude[i].values)] for i in range(len(ds9.id_dim))]
list_coords_ds7 = [ [float(ds7.longitude[i].values), float(ds7.latitude[i].values)] for i in range(len(ds7.id_dim))]


uniq7, xx7 ,yy7, ind7 = unique_coords(list_coords_ds7)
uniq9, xx9 ,yy9, ind9 = unique_coords(list_coords_ds9)

print(f" Size of co7: {len(ds7.id_dim)}. Number of unique pts: {len(uniq7)}")
print(f" Size of co9: {len(ds9.id_dim)}. Number of unique pts: {len(uniq9)}")


########################
# Find closest neighbours in co9 for each co7 member

def find_closest_ind( xpt,ypt, xarr,yarr, verbose=False ):
    """ find the index of the closest point in the array of coordinates """
    ind = np.argmin( np.square(xarr - xpt).values + np.square(yarr - ypt).values )
    if verbose: print(f"ind:{ind} xpt:{xpt.values} xarr[i]:{xarr[ind].values}.  ypt:{ypt.values} yarr[i]:{yarr[ind].values}")
    return ind

## Check how np.unique works
#np.unique( [1,2,3,4,4,3,2], return_index=True, return_counts=True)

# Accumulate the indices from co9 that match co7 points
II_co9_from_co7 = []
for i in range(len(ds7.id_dim)):
    II_co9_from_co7.append(find_closest_ind( ds7.longitude[i], ds7.latitude[i], ds9.longitude, ds9.latitude, verbose=True) )
    

print(f"The indices from CO9 that match CO7: {II_co9_from_co7}")

########################

# Define some short hand and plot to visually confirm

ds9['lon'] = np.mod(ds9.longitude+180, 360)-180
ds7['lon'] = np.mod(ds7.longitude+180, 360)-180

xx7 = ds7['lon']
yy7 = ds7['latitude']

xx9 = ds9['lon'][II_co9_from_co7]
yy9 = ds9['latitude'][II_co9_from_co7]


########################


fig, [ax1,ax2] = plt.subplots(1, 2, figsize=(10,8))

ax1.plot(xx7, yy7, 'bx')
ax1.plot(xx9, yy9, 'r+')
for i in range(len(ds7.id_dim)):
 ax1.text(xx7[i], yy7[i], i, color='b')
ax1.set_title('unique indices from co7')
ax1.legend(['co7','co9'])

ax2.plot(xx7, yy7, 'bx')
ax2.plot(xx9, yy9, 'r+')
for i in range(len(ds7.id_dim)):
 ax2.text(xx9[i], yy9[i], II_co9_from_co7[i], color='r')
ax2.set_title('unique indices from co9')


########################

fig, [ax1,ax2] = plt.subplots(1, 2, figsize=(10,8))

ax1.plot(np.sort(xx7), np.sort(xx9), '+')
ax1.plot([-10,10],[-10,10])

ax2.plot(np.sort(yy7), np.sort(yy9), '+')
ax2.plot([48,62],[48,62])

########################
# Save files with reduced set of matching ports

ds9.isel(id_dim=II_co9_from_co7).to_netcdf(fn_co9.replace('P1.5c','P1.5c_matchCO7'))

