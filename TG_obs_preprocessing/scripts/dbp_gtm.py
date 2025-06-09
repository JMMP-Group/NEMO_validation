#For GlobalTideModel specific routines, e.g. plotting. Probably easily transferrable to other
#projects.

import numpy as np
from netCDF4 import Dataset as ds
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cmocean.cm as cmo
import scipy.ndimage.filters as spf
import cartopy.crs as ccrs
import cartopy as cartopy
import matplotlib.colors as colors
import dbp_general as dbg
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import scipy.interpolate as intp
import os
import time
#from global_land_mask import globe
import multiprocessing as mp

def ensemble_correlation(ens1, ens2):
    
    ne, nr, nc = ens1.shape
    ens1 = ens1 - np.nanmean(ens1, axis=0)
    ens2 = ens2 - np.nanmean(ens2, axis=0)
    
    c = np.zeros((nr,nc))*np.nan
    
    for ii in range(0,nr):
        for jj in range(0,nc):
            c[ii,jj] = np.corrcoef(ens1[:,ii,jj], ens2[:,ii,jj])[0,1]
    
    return c

def create_land_mask_par(lon, lat, n_proc = 24):
    '''
    # Runs the create_land_mask() routine in parallel using 
    # multiprocessing.starmap. 
    '''
    
    lonF = lon.flatten()
    latF = lat.flatten()
    print('Starting search..')
    
    l180 = lonF>180
    lonF[l180] = lonF[l180] - 360
    
    n_pts = len(lonF)
    landmaskF = np.zeros(n_pts)
    
    print('Starting parallel opts..')
    pool = mp.Pool(n_proc)
    landmaskF = pool.starmap(globe.is_land, [(latF[ii], lonF[ii]) for ii in range(0,n_pts)])
    pool.close()
    pool.join()
    print('Pool closed...')
    
    print('reshaping..')
    landmask = np.reshape(landmaskF, lon.shape)
    
    return landmask

def create_land_mask(lon, lat):
    '''
    # Creates a 2-dimensional land mask from 2-dimension geolocation data.
    # Uses global-land-mask functionality but basemap can also be used too
    # (bm.is_land). Can take a long time for larger datasets. In
    # this case, use the parallel routine create_land_mask_par.
    '''
    
    lonF = lon.flatten()
    latF = lat.flatten()
    print('Starting search..')
    landmaskF = [globe.is_land(latF[ii], lonF[ii]) for ii in range(0,len(lonF))]
    print('reshaping..')
    landmask = np.reshape(landmaskF, lon.shape)
    
    return landmask

def GTM_file_create_obs(fn, n_obs, c):
    '''
    # Creates a new GTM obs file, overwriting an old one if it exists. The base 
    # file contains just the locations dimension loc. 
    '''
    # Delete existing file (if exists)
    try:
        os.system('rm -f '+ fn)
        print('Overwriting existing file.')
    except:
        print('Creating new file')
    
    # Open file for writing
    ncid = ds(fn,'w')
    
    # Get data size and create dimensions
    dim_obs = ncid.createDimension('locs', n_obs)
    
    ncid.title = 'Observations used for GTM DA for ' + c + ' constituent.' 
    
    ncid.history = 'Created ' + time.ctime(time.time())
    ncid.close()
    
    return 

def GTM_file_append_obs(fn, var, varname, varunits = ''):
    '''
    # Adds a variable to an existing GTM obs file. 
    '''
    
    # Open file for appending
    ncid = ds(fn, 'a')
    varnc = ncid.createVariable(varname, np.float32, ('locs'))
    varnc[:] = var
    varnc.units = varunits
    ncid.close()
    
    return 

def GTM_file_create_ungridded(fn, lon, lat, c, datatype='z'):
    '''
    # Creates a new GTM file, overwriting an old one if it exists. The base 
    # file contains just lon and lat arrays (so that dimensions can be 
    # established) and has two datatype options: z and uv. These options 
    # just control the description options in the file. c is the constituent
    '''
    # Delete existing file (if exists)
    try:
        os.system('rm -f '+ fn)
        print('Overwriting existing file.')
    except:
        print('Creating new file')
    
    # Open file for writing
    ncid = ds(fn,'w')
    
    # Get data size and create dimensions
    nr,nc = lon.shape
    dim_y = ncid.createDimension('y', nr)
    dim_x = ncid.createDimension('x', nc)
    
    # Create lon/lat variables
    var_lat = ncid.createVariable('latitude' , np.float32, ('y','x'))
    var_lon = ncid.createVariable('longitude', np.float32, ('y','x'))
    var_lat[:,:] = lat
    var_lon[:,:] = lon
    
    if datatype=='z':
        ncid.title = 'NOCGTM (ungridded) v1.0 | Heights for ' + c + ' constituent.' 
    
    ncid.history = 'Created ' + time.ctime(time.time())
    ncid.close()
    
    return print('Written')

def GTM_file_append_ungridded(fn, auv, varname, varunits = ''):
    '''
    # Adds a variable to an existing GTM .nc file. Needs to be orientated the
    # same way as for the original file. 
    '''
    
    # Open file for appending
    ncid = ds(fn, 'a')
    var = ncid.createVariable(varname, np.float32, ('y','x'))
    var[:,:] = auv
    var.units = varunits
    ncid.close()
    
    return print('Written')

def GTM_file_create_gridded(fn, lon, lat, c, datatype='z'):
    '''
    # Creates a new GTM file, overwriting an old one if it exists. The base 
    # file contains just lon and lat arrays (so that dimensions can be 
    # established) and has two datatype options: z and uv. These options 
    # just control the description options in the file. c is the constituent
    '''
    # Delete existing file (if exists)
    try:
        os.system('rm -f '+ fn)
        print('Overwriting existing file.')
    except:
        print('Creating new file')
    
    # Open file for writing
    ncid = ds(fn,'w')
    
    # Get data size and create dimensions
    nr = len(lon)
    nc = len(lat)
    dim_x = ncid.createDimension('latitude', nc)
    dim_y = ncid.createDimension('longitude', nr)
    
    # Create lon/lat variables
    var_lat = ncid.createVariable('latitude' , np.float32, 'latitude')
    var_lon = ncid.createVariable('longitude', np.float32, 'longitude')
    var_lat[:] = lat
    var_lon[:] = lon
    
    if datatype=='z':
        ncid.title = 'NOCGTM v1.0 | Heights for ' + c + ' constituent.' 
    else:
        ncid.title = 'NOCGTM v1.0 | Currents for ' + c + ' constituent.'
    
    ncid.history = 'Created ' + time.ctime(time.time())
    ncid.close()
    
    return print('Written')

def GTM_file_append_gridded(fn, auv, varname, varunits = ''):
    '''
    # Adds a variable to an existing GTM .nc file. Needs to be orientated the
    # same way as for the original file. 
    '''
    
    # Open file for appending
    ncid = ds(fn, 'a')
    var = ncid.createVariable(varname, np.float32, ('latitude','longitude'))
    var[:,:] = auv
    var.units = varunits
    ncid.close()
    
    return print('Written')

def GTM_grid_variable(x, y, v, x2, y2, method='nearest'):
    '''
    # Interpolation routine to grid NEMO ORCA grid data to a regular 1/12 deg
    # grid. 
    '''
    #shift longitudes
    x[x>180] = x[x>180] - 360
    
    #copy arrays
    x = np.array( x ); y = np.array( y ); v = np.array( v )
    #flatten arrays
    x = x.flatten(); y = y.flatten(); v = v.flatten()
    #apply landmask (remove land points)
    mask = np.isnan(v)
    x = x[~mask]; y = y[~mask]; v = v[~mask]
    print('Done processing')
    
    xgrid, ygrid = np.meshgrid(x2,y2)
    points = (x, y)
    vint = intp.griddata(points , v, (xgrid, ygrid), method=method)
    
    return vint

def convert_nc_to_textfile(fn_in, dn_out, const, mask):
    
    ncid = ds(fn_in)
    
    lat = ncid.variables['latitude'][:,:].flatten()
    lon = ncid.variables['longitude'][:,:].flatten()
    mask = mask.flatten()
    
    for cc in const:
        a = ncid.variables['M2_a'][:,:].flatten()
        g = ncid.variables['M2_g'][:,:].flatten()
        ua = ncid.variables['M2_ua'][:,:].flatten()
        ug = ncid.variables['M2_ug'][:,:].flatten()
        va = ncid.variables['M2_va'][:,:].flatten()
        vg = ncid.variables['M2_vg'][:,:].flatten()
        
        a = a[~mask]
        g = g[~mask]
        ua = ua[~mask]
        ug = ug[~mask]
        va = va[~mask]
        vg = vg[~mask]
        
        L = len(a)
        
        fn_tmp = 'GTM_' + cc + '.dat'
        fn_out = dn_out + '/' + fn_tmp
        print(fn_out)
        fid = open(fn_out, 'w')
        for ii in range(0,L,2):
            tmp1 = str(lat[ii]) + ' ' + str(lon[ii]) + ' '
            tmp2 = str(a[ii]) + ' ' + str(g[ii]) + ' '
            tmp3 = str(ua[ii]) + ' ' + str(ug[ii]) + ' '
            tmp4 = str(va[ii]) + ' ' + str(vg[ii]) + ' '
            line = '{0:9.3f} {1:9.3f} {2:9.3f} {3:9.3f} {4:9.3f} {5:9.3f} {6:9.3f} {7:9.3f}'
            line.format(lat[ii], lon[ii], a[ii], g[ii], ua[ii], ug[ii], va[ii], vg[ii])
            print(line)
            fid.write(line)
        fid.close()
    
    return print('done')

def find_local_minima2(X):
    '''
    # Finds local minima in a 2D array. Output is a pair of indices (row,col)
    #
 IN # X     :: Input 2-dimensional array.
    #
OUT # 
    '''
    
    neighborhood_size = 5
    threshold = 0.1
    
    X_max = filters.maximum_filter(X, neighborhood_size)
    maxima = (X == X_max)
    X_min = filters.minimum_filter(X, neighborhood_size)
    minima = (X == X_min)
    diff = ((X_max - X_min) > threshold)
    minima[diff == 0] = 0

    labeled, num_objects = ndimage.label(minima)
    slices = ndimage.find_objects(labeled)
    x, y = [], []
    for dy,dx in slices:
        x_center = (dx.start + dx.stop - 1)/2
        x.append(x_center)
        y_center = (dy.start + dy.stop - 1)/2    
        y.append(y_center)
        
    return x, y

def find_local_maxima2(X, neighborhood_size = 5, threshold = 75):
    '''
    # Finds local maxima in a 2D array. Output is a pair of indices (row,col)
    #
 IN # X     :: Input 2-dimensional array.
    #
OUT # 
    '''
    
    X_max = filters.maximum_filter(X, neighborhood_size)
    maxima = (X == X_max)
    X_min = filters.minimum_filter(X, neighborhood_size)
    diff = ((X_max - X_min) > threshold)
    maxima[diff == 0] = 0

    labeled, num_objects = ndimage.label(maxima)
    slices = ndimage.find_objects(labeled)
    x, y = [], []
    for dy,dx in slices:
        x_center = (dx.start + dx.stop - 1)/2
        x.append(x_center)
        y_center = (dy.start + dy.stop - 1)/2    
        y.append(y_center)
        
    return x, y


def identify_amphidromes(lon, lat, g, neighborhood_size = 5, threshold = 75):
    
    lon = lon[5:-5,:]
    lat = lat[5:-5,:]
    g = g[5:-5,:]
    
    nr, nc = np.shape(lon)
    
    #shift_lon = np.zeros((nr,nc))
    #shift_lat = np.zeros((nr,nc))
    
    shift_lon = 0.5*(lon[:,:-1] + lon[:,1:])
    shift_lat = 0.5*(lon[:-1,:] + lon[1:,:])
    
    dx = 1
    dy = 1
    
    Ag_x = np.zeros((nr,nc))
    Ag_y = np.zeros((nr,nc))
    
    Ag_x[:,0:-1] = dbg.compare_angles(g[:,0:-1], g[:,1:])/dx
    Ag_y[0:-1,:] = dbg.compare_angles(g[0:-1,:], g[1:,:])/dy
    
    Ag_s = np.sqrt(Ag_x**2 + Ag_y**2)
    
    maxi, maxj = find_local_maxima2(Ag_s, neighborhood_size = 5, 
                                    threshold = 75)
    #y_maxi, y_maxj = find_local_maxima2(Ag_y)
    
    maxi = np.array(maxi).astype(int)
    maxj = np.array(maxj).astype(int)
    #y_maxi = np.array(y_maxi).astype(int)
    #y_maxj = np.array(y_maxj).astype(int)
    
    amphi_lon = lon[maxj,maxi]
    amphi_lat = lat[maxj,maxi]
    
    return amphi_lon, amphi_lat


def gen_regular_grid_file(filename, X, Y, mask):
    
    ncid = ds(filename,'w')
    
    dim_x = ncid.createDimension('x',np.shape(X)[1])
    dim_Y = ncid.createDimension('y',np.shape(X)[0])
    
    nc_lat = ncid.createVariable('lat','f4',('y','x'))
    nc_lon = ncid.createVariable('lon','f4',('y','x'))
    nc_mask = ncid.createVariable('mask','f4',('y','x'))
    
    nc_lon[:,:] = X[:,:]
    nc_lat[:,:] = Y[:,:]
    nc_mask[:,:] = mask[:,:]
    
    ncid.close()
    
    return

def smooth_bathy(file_in, file_out, sigma=1, applymask=0, cut=0, mindep=10):
    #Reads in and smooths a bathymetry file.
    
    ncin = ds(file_in)
    
    bathy = ncin.variables['Bathymetry'][:,:]
    lon   = ncin.variables['nav_lon'][:,:]
    lat   = ncin.variables['nav_lat'][:,:]
    
    ds.close(ncin)
    
    bathy = spf.gaussian_filter(bathy,sigma)
    
    if cut==1:
        logic1 = bathy<mindep
        logic2 = bathy==0
        logic3 = logic1*(~logic2)
        bathy[logic3] = mindep
    ncout = ds(file_out,'w')
   
    n_row, n_col = np.shape(bathy)
    
    dim_x = ncout.createDimension('x',n_col)
    dim_y = ncout.createDimension('y',n_row)
    
    var_bath = ncout.createVariable('Bathymetry','f4',('y','x'))
    var_lat = ncout.createVariable('nav_lat','f4',('y','x'))
    var_lon = ncout.createVariable('nav_lon','f4',('y','x'))
    
    var_bath[:,:] = bathy
    var_lat[:,:] = lat
    var_lon[:,:] = lon
    
    ncout.close()
    
    return bathy

def gen_constant_bath_file(depth, in_file, out_file):
    #Generates a constant bathymetry file suited for use with NEMO. Depth should
    #be positive. in_file is a NEMO coordinates file. out_file is the output
    #file (e.g. bathy_meter.nc)
    
    nc_in = ds(in_file,'r')
    nav_lon = nc_in.variables['nav_lon'][:,:]
    nav_lat = nc_in.variables['nav_lat'][:,:]
    n_row,n_col = np.shape(nav_lon)
    
    nc_in.close()
    
    depth = np.zeros((n_row,n_col)) + depth
    
    nc_out = ds(out_file,'w')
    
    dim_x = nc_out.createDimension('x',n_col)
    dim_y = nc_out.createDimension('y',n_row)
    
    var_bath = nc_out.createVariable('Bathymetry','f4',('y','x'))
    var_lat = nc_out.createVariable('nav_lat','f4',('y','x'))
    var_lon = nc_out.createVariable('nav_lon','f4',('y','x'))
    
    print(nav_lat[1,1])
    
    var_bath[:,:] = depth
    var_lat[:,:] = nav_lat
    var_lon[:,:] = nav_lon
    
    nc_out.close()
    
    return depth

def set_cartopy_domain(domain):
    '''
    # Set domain and projection for cartopy plots. Function output can be used
    # directly with plotting functions if cartopy is loaded.
    '''
    
    labels = 0
    
    if domain=='global_miller': #Global Miller Projection
        ax = plt.axes(projection=ccrs.Miller(central_longitude=0))
        ax.set_global()
    elif domain == 'global_robinson': #Global Robinson Projection
        ax = plt.axes(projection=ccrs.Robinson(central_longitude=0))
        ax.set_global()
    elif domain == 'global_mercator': #Global Robinson Projection
        ax = plt.axes(projection=ccrs.Mercator(central_longitude=0))
        ax.set_global()
        labels = 1
        
    print("Set domain to " + domain)
    
    return ax, labels

def tide_map_gtm(lons, lats, amp, phase, ampmax=-999, title='', 
                 domain = 'global_miller', gridlines=0):
    '''
    # For creating a global tide map with amplitude filled contours and phase
    # contours. 
    #
 IN # lons, lats :: 2D longitude/latitude grids
 IN # amp, phase :: 2D amplitude/phase grids
 IN # ampmax     :: maximum amplitude for colorbar
 IN # gridlines  :: whether or not to display lat/lon gridlines
    #
OUT # None
    '''
    
    plt.figure(figsize=(10,7))
    n_levels=50
    
    if ampmax == -999:
        ampmax = np.nanmax(amp)

    ax, labels = set_cartopy_domain(domain)

    levels = np.linspace(0,ampmax,n_levels)
    norm = colors.BoundaryNorm(levels, ncolors=256)
    
    mesh = ax.pcolormesh(lons, lats, amp, cmap='cmo.amp',norm=norm, 
                         transform=ccrs.PlateCarree())
    plt.colorbar(mesh,orientation='vertical', shrink=0.8, extend='max')
    mesh = ax.contour(lons[2:-2,:-2], lats[2:-2,:-2], phase[2:-2,:-2],
                      transform=ccrs.PlateCarree(), colors='k', 
                      linewidths=0.75, levels = 10)
    if gridlines:
        if labels:
            ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
            gl.xlabels_top = False
            gl.ylabels_right = False
        else:
            ax.gridlines()
    ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.75)
    plt.title(title)
    plt.show()
    
    return

def geo_scatter(lons, lats, C = 'b', s = 10, marker = ['o'], title='', 
                domain='global_miller', gridlines=True, cmap='cmo.amp', 
                landcolor = [0.9, 0.9, 0.9], clims = 'auto'):
    '''
    # Plots scattered points on a geographic cartopy plot. Can either plot a
    # single pair of datasets or multiple by supplying specified inputs as
    # tuples
    #
 IN # lons, lats   :: 1D longitude/latitude arrays or multiple 1D arrays inside
    #                 a tuple
 IN # C            :: Colors for scattered points. Can be vector or single
    #                 color string or RGB. Or multiple inside a tuple. If fewer
    #                 colors than datasets are supplied, they will loop
 IN # s, marker    :: Marker size and type (type is a list if multiple 
    #                 datasets)
 IN # domain       :: Choice of plotting domain. See set_cartopy_domain
 IN # gridlines    :: If true, plot lat/lon grid lines
 IN # cmap         :: colormap string (get_cmap) or specified values
    #
OUT # fig, ax      :: Matplotlib figure/axis objects for plot.
    '''
    
    if type(lons) is tuple:
        n_ds = len(lons)
        n_color = len(C)
        n_marker = len(marker)
    else:
        n_ds = 1
        C = (C)
        n_color = 1
        n_marker = 1
    
    if type(cmap) == str:
        cmap = cm.get_cmap(cmap)
    
    f = plt.figure(figsize = (10,7))
    ax, labels = set_cartopy_domain(domain)
    
    for ii in range(0,n_ds):
        if n_ds == 1:
            print(C[ii%n_color])
            sca = ax.scatter(lons, lats, c=C, s=s, zorder=10, 
                             marker=marker[ii%n_marker],
                             transform=ccrs.PlateCarree(), cmap=cmap)
        else:
            sca = ax.scatter(lons[ii], lats[ii], c=C[ii%n_color], s=s, 
                             zorder=10, marker=marker[ii%n_marker],
                             transform=ccrs.PlateCarree(), cmap=cmap)
            
    plt.colorbar(sca, orientation='vertical', shrink=0.8, extend='max')
        
    ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
    ax.add_feature(cartopy.feature.LAND,color=landcolor)
    if gridlines:
        if labels:
            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
            gl.xlabels_top = False
            gl.ylabels_right = False
        else:
            ax.gridlines()
    #ax.colorbar(extend='max')
    plt.title(title)
    plt.show()
    
    return f, ax

def contourf_gtm(lons, lats, amp, ampmax=-999, title='', 
                 domain = 'global_miller', gridlines=0):
    '''
    # For creating a global tide map with amplitude filled contours and phase
    # contours. 
    #
 IN # lons, lats :: 2D longitude/latitude grids
 IN # amp, phase :: 2D amplitude/phase grids
 IN # ampmax     :: maximum amplitude for colorbar
 IN # gridlines  :: whether or not to display lat/lon gridlines
    #
OUT # None
    '''
    
    plt.figure(figsize=(10,7))
    n_levels=50
    
    if ampmax == -999:
        ampmax = np.nanmax(amp)

    ax, labels = set_cartopy_domain(domain)

    levels = np.linspace(0,ampmax,n_levels)
    norm = colors.BoundaryNorm(levels, ncolors=256)
    
    mesh = ax.pcolormesh(lons, lats, amp, cmap='cmo.amp', norm=norm, 
                         transform=ccrs.PlateCarree())
    plt.colorbar(mesh, orientation='vertical', shrink=0.8, extend='max')

    if gridlines:
        if labels:
            ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
            gl.xlabels_top = False
            gl.ylabels_right = False
        else:
            ax.gridlines()
    ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.75)
    plt.title(title)
    plt.show()
    
    return

def discrete_pd2(A,dx=1,dy=1):
    #Quick estimation of horizontal derivatives in x and y.
    
    nr, nc = np.shape(A)
    Ag_x = np.zeros((nr,nc))
    Ag_y = np.zeros((nr,nc))
    
    Ag_x[:,0:-1] = np.abs(A[:,0:-1] - A[:,1:])/dx
    Ag_y[0:-1,:] = np.abs(A[0:-1,:] - A[1:,:])/dy
    
    return Ag_x, Ag_y

def discrete_pa2(A,dx=1,dy=1):
    #Quick estimation of horizontal derivatives in x and y.
    
    nr, nc = np.shape(A)
    Ag_x = np.zeros((nr,nc))
    Ag_y = np.zeros((nr,nc))
    
    #Ag_x[:,0:-1] = np.abs(A[:,0:-1] - A[:,1:])/dx
    #Ag_y[0:-1,:] = np.abs(A[0:-1,:] - A[1:,:])/dy
    
    Ag_x[:,0:-1] = dbg.compare_angles(A[:,0:-1], A[:,1:])/dx
    Ag_y[0:-1,:] = dbg.compare_angles(A[0:-1,:], A[1:,:])/dy
    
    return Ag_x, Ag_y
