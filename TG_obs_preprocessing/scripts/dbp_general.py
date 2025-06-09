"""
Created on Tue Feb 12 14:25:25 2019

@author: dbyrne
A set of general-use functions for general and mathsy things in general
"""

import numpy as np
from netCDF4 import Dataset as ds
import datetime as datetime
import scipy.interpolate as spi
from scipy.optimize import curve_fit

def cov_partial(var1, var2):
    '''
    # Calculates a part of covariance matrix consisting of over cross-variable
    # pairs. I.E how does each variables in var1 covary with each variable in
    # var2? Input should be of form n_var*n_obs, with variables on rows. Var2
    # will be on the rows of the resulting covariance matrix.
    '''
    nv1, no1 = var1.shape
    nv2, no2 = var2.shape
    
    v1mean = np.nanmean(var1, axis=1)
    v2mean = np.nanmean(var2, axis=1)
    
    
    C = np.zeros()
    
    return C

def gauss(x, sigma):
        p = [sigma]
        return np.exp(-((x)/p)**2)
    
def surface1(xy, a, b, c, d):
    x, y = xy
    out = a + b*x + c*y + d*x*y
    return np.ravel(out)

def fit_gauss(x, y):
    popt, pcov = curve_fit(gauss, x, y)
    return popt, pcov

def ext_var_ind_box(rbounds, cbounds, var_rr, var_cc):
    '''
    # For extracting elements of variable vector that are within a index box
    # defined by rbounds (rows) and cbounds (cols). Output is an array of 
    # indices.
    '''
    var_rr = np.array(var_rr)
    var_cc = np.array(var_cc)
    logr1 = var_rr >= rbounds[0]
    logr2 = var_rr <= rbounds[1]
    logc1 = var_cc >= cbounds[0]
    logc2 = var_cc <= cbounds[1]
    logr = logr1*logr2
    logc = logc1*logc2
    box_ind = np.where(logr*logc == 1)    
    return box_ind

def sub2ind(array_shape, rows, cols):
    '''
    # Converts 2D array index to 1D location assuming an array has been 
    # flattened along the first dimension (row-major)
    '''
    array_shape = np.array(array_shape)
    rows = np.array(rows)
    cols = np.array(cols)
    return rows*array_shape[1] + cols

def angular_distance(lon1, lon2, lat1, lat2):
    '''
    Angular distance between two geographical points
    '''
    lon1 = np.radians(lon1)
    lon2 = np.radians(lon2)
    lat1 = np.radians(lat1)
    lat2 = np.radians(lat2)
    
    s1 = np.sin(lon1)*np.sin(lon2) + np.cos(lon1)*np.cos(lon2)*np.cos(lat1-lat2)
    alpha = np.arccos(s1)
    
    return np.degrees(alpha)

def r_squared_lin(x, y, fit):
    '''For calculating r-squared of a linear fit. Fit should be a python polyfit
    object'''
    
    fity = fit(x)
    diff = (y - fity)**2
    ybar = np.nanmean(y)
    ymybar = (y - ybar)**2
    
    SStot = np.nansum(ymybar)
    SSres = np.nansum(diff)
    
    R2 = 1 - SSres/SStot
    
    return R2

def dist_matrix(x, y):
    '''
    # Calculates a distance matrix on a rectangular grid (pythagoras).
    # Matrix will contain distances between every pair of points. 
    #
 IN # x1, y1 :: First set of locations.
 IN # x2, y2 :: Second set of locations.
    '''
    x = np.array(x)
    y = np.array(y)

    n_data = len(x)
    
   #Define distance array D.
    D = np.zeros((n_data, n_data))
    
    #Loop over every pair of points and calculate distance using pythagoras
    for rr in range(0, n_data):
        for cc in range(0, n_data):
            D[rr,cc] = dist_pythagoras(x[rr], y[rr], x[cc], y[cc])
    
    return D

def dist_pythagoras(x1, y1, x2, y2):
    '''
    Calculates pythagorean distance between two points (x1, y1) and (x2, y2)
    '''
    xdist = (x1 - x2)**2
    ydist = (y1 - y2)**2
    xysum = xdist + ydist
    dist = np.sqrt(xysum)
    
    return dist

def dist_matrix_haversine(lon1, lat1, lon2 = None, lat2 = None, degrees = True, 
                          R = 6371.00718, alg = 'part'):
    '''
    # Calculates a geographic autodistance matrix using the haversine function.
    # Matrix will contain distances between every pair of points. This is 
    # slow at the moment and definitely needs optimising. The second set of
    # locations will end up as the rows of the distance matrix.
    #
 IN # lon:: Vector of longitudes.
 IN # lon:: Vector of Latitudes.
 IN # degrees    :: If true then will be converted  to radians.
 IN # R          :: Radius of the earth in desired output units (default km).
    '''
    L1 = len(lon1)
    L2 = len(lon2)
    
    if alg == 'sym':
        X1, X2 = np.meshgrid(lon1, lon1)
        Y1, Y2 = np.meshgrid(lat1, lat1)
    elif alg == 'part':
        X1 = np.tile(lon1,[L2,1])
        Y1 = np.tile(lat1,[L2,1])
        X2 = np.transpose(np.tile(lon2,[L1,1]))
        Y2 = np.transpose(np.tile(lat2,[L1,1]))

    D = dist_haversine(X1,Y1,X2,Y2, degrees=degrees)
    
    return D


def dist_matrix_sphere(lat,lon, R=6371.007176 ):
    """
    Compute a distance matrix of the coordinates using a spherical metric.
    :param coordinate_array: numpy.ndarray with shape (n,2); latitude is in 
    1st col, longitude in 2nd. :returns distance_mat: numpy.ndarray with shape 
    (n, n) containing distance in km between coords.
    """

    latitudes = lat
    longitudes = lon

    # Convert latitude and longitude to spherical coordinates in radians.
    degrees_to_radians = np.pi/180.0
    phi_values = (90.0 - latitudes)*degrees_to_radians
    theta_values = longitudes*degrees_to_radians

    # Expand phi_values and theta_values into grids
    theta_1, theta_2 = np.meshgrid(theta_values, theta_values)
    theta_diff_mat = theta_1 - theta_2

    phi_1, phi_2 = np.meshgrid(phi_values, phi_values)

    # Compute spherical distance from spherical coordinates
    angle = (np.sin(phi_1) * np.sin(phi_2) * np.cos(theta_diff_mat) + 
           np.cos(phi_1) * np.cos(phi_2))
    arc = np.arccos(angle)

    # Multiply by earth's radius to obtain distance in km
    return arc * R

def dist_haversine(lon1, lat1, lon2, lat2, degrees=True, R = 6371.007176 ):
    '''
    # Estimation of geographical distance using the Haversine function.
    # Input can be single values or 1D arrays/lists of locations. This
    # does NOT create a distance matrix but outputs another 1D array.
    # This works for either location vectors of equal length OR a single loc
    # and an arbitrary length location vector.
    #
    # lon1, lat1 :: Start location(s).
    # lon2, lat2 :: End location(s).
    # degrees    :: If true then will be converted  to radians.
    # R          :: Radius of the earth in desired output units (default km).
    
    #Convert to numpy array if necessary (if single number, it should be fine)
    '''
    if type(lon1) == list:
        lat1 = np.array(lat1)
        lon1 = np.array(lon1)
        lat2 = np.array(lat2)
        lon2 = np.array(lon2)
    #Convert to radians if necessary
    if degrees:
        lat1 = np.deg2rad(lat1)                     
        lon1 = np.deg2rad(lon1)     
        lat2 = np.deg2rad(lat2)                       
        lon2 = np.deg2rad(lon2)  
    
    #Latitude and longitude differences
    dlat = (lat2-lat1)/2
    dlon = (lon2-lon1)/2
    
    #Haversine function.
    d = np.sin(dlat)**2 + np.cos(lat1)*np.cos(lat2) * np.sin(dlon)**2
    d = 2 * R * np.arcsin(np.sqrt(d))
    
    return d 

def geo_nn_indices(lon1, lat1, lon2, lat2, mask=None, 
                   degrees=True, R = 1):
    '''
    # Determines indices for a NEAREST NEIGHBOUR interpolation from 
    # (lon1, lat1) to (lon2, lat2). Interpolation can be done by extracting 
    # the output indices from the input array data. This uses the haversine
    # function to estimate spherical distances (rather than euclidean). Output
    # is nr_ii (row index) and nr_jj (column index).
    #
    # lon1, lat1 :: Location data 1 (interpolating from these points). Should
    #               be 2-dimensional.
    # lon2, lat2 :: Location data 2 (interpolating to these points), should be
    #               1-dimensional.
    # mask       :: Exclude True points from location data 1 (e.g. land points)
    # degrees    :: If true then will be converted  to radians.
    # R          :: Radius of the earth in desired output units (default km).
    '''
    
    #Ensure input arrays are numpy arrays
    lon1 = np.array(lon1); lat1 = np.array(lat1)
    lon2 = np.array(lon2); lat2 = np.array(lat2)
    nr, nc = np.shape(lon1)
    #Define index arrays for later indexing
    X, Y = np.meshgrid(np.arange(0,nc), np.arange(0,nr))
    
    #Reduce input arrays using mask (if required), otherwise flatten
    if type(mask) is None:
        lon1 = lon1.flatten(); lat1 = lat1.flatten()
        X = X.flatten(); Y = Y.flatten()
    else:
        lon1 = reduce_array_using_mask(lon1, mask)
        lat1 = reduce_array_using_mask(lat1, mask)
        X = reduce_array_using_mask(X, mask)
        Y = reduce_array_using_mask(Y, mask)

    #Initialisation of output arrays
    nr_ii = np.zeros(len(lon2))
    nr_jj = np.zeros(len(lat2))
    
    #D = dist_matrix_haversine(lon1, lat1, lon2, lat2)
    #print('x')
    #minii = np.argmin(D, axis=1)
    #nr_ii = Y[minii]
    #nr_jj = X[minii]
    
    for ii in range(0,len(lon2)):
        dist = dist_haversine(lat2[ii], lon2[ii],
                                  lat1, lon1, degrees=degrees, R=R)        
        mini = np.argmin(dist)
        nr_ii[ii] = Y[mini]
        nr_jj[ii] = X[mini]
    
    return nr_ii.astype(int), nr_jj.astype(int)

def rebuild_array_using_mask(A_red, mask):
    '''
    # Rebuilds an array that was reduced using reduce_array_using_mask. A_red
    # is a 1-dimensional array and mask is the 2-dimensional array used to
    # reduce the original array. Output is 2-dimensional.
    '''
    return A

def reduce_array_using_mask(A, mask):
    '''
    # Removes points from an array A based on values in mask. Input arrays A
    # and mask are 2-dimensional and output is a 1-dimensional reduced array.
    # Reduction is done vertically (along columns first. Hence 'F'). It is the True values
    # in the mask that are removed from A.
    '''
    A = np.array(A)
    A = A.flatten('F')
    mask = np.array(mask)
    mask = mask.flatten('F')
    A_red = A[ np.where(mask==0) ]
    return A_red

def date_parser(year, month, day, hour=0,minute=0,second=0):
    year, month, day, hour = map(int, (year, month, day, hour))
    return datetime.datetime(year, month, day, hour, minute, second)

def determine_ranges(max_val, max_ii, min_val, min_ii):
    #Calculates difference between successive maxima and minima of a time 
    #series already determined using find_optima.
    
    #if max_ii[0] < min_ii[0]:
    ranges = max_val - min_val
    
    return ranges

def apply_fill(A, fillvalues):
    # Applies a list of fill values to an array A. Values specific in
    # fillvalues list will be change to np.nan.
    
    if type(fillvalues) != list:
        fillvalues = [fillvalues]
    
    A = np.array(A)
    for ii in range(0,len(fillvalues)):
        A[A == fillvalues[ii]] = np.nan
    
    return A

def apply_mask(A,mask):
    # Applies 2D mask to a 2D array A. Replaces locations where mask = 1 with 
    # nan 
    A = np.array(A)
    A[mask] = np.nan
    return A

def compare_angles(a1,a2,degrees=True):
    #Compares the difference between two angles. e.g. it is 2 degrees between
    #359 and 1 degree. If degrees = False then will treat angles as radians.
    
    if not degrees:
        a1 = np.degrees(a1)
        a2 = np.degrees(a2)
        
    diff = 180 - np.abs(np.abs(a1-a2)-180)
    
    if not degrees:
        a1 = np.radians(a1)
        a2 = np.radians(a2)
    
    return diff

def cart2polar(x, y, degrees = True):
    '''
    # Conversion of cartesian to polar coordinate system
    # Output theta is in radians
    '''
    r     = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    if degrees:
        theta = np.rad2deg(theta)
    return(r, theta)

def polar2cart(r, theta, degrees=True):
    '''
    # Conversion of polar to cartesian coordinate system
    # Input theta must be in radians
    '''
    if degrees:
        theta = np.deg2rad(theta)
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return(x, y)

def convert_harmonics(z1,z2):
    #Converts harmonics from complex space to amp/phase space. Output in degrees.
    #Will work for scalars or arrays (elementwise). 
    
    n_rows,n_cols = np.shape(z1)
    
    amp, phase = cart2polar(z1,z2)
    phase = np.rad2deg(phase)
    
    for ii in range(0,n_rows):
        for jj in range(0,n_cols):
            if phase[ii,jj] < 0:
                phase[ii,jj] = phase[ii,jj] + 360
    
    return amp, phase

def rmsd(x1,x2):
    #Calculated rmsd of 2 given vectors of values x1 and x2. One is likely obs
    #and the other model output.
    
    diff = x1 - x2
    diffs = diff**2
    msd = np.nanmean(diffs)
    rmsd = np.sqrt(msd)
    
    return rmsd

def rmsd2(x1,x2):
    #calculates RMSD for each row of arrays x1 and x2
    
    n_rows = np.shape(x1)[0]
    rmsd_out = np.empty(n_rows)
    
    for ii in range(0,n_rows):
        rmsd_out[ii] = rmsd(x1[ii,:],x2[ii,:])
    
    return rmsd_out

def md(x1,x2):
    #Calculates the mean differences between two time series x1 and x2
    
    diff = x1 - x2
    md = np.nanmean(diff)
    
    return md

def md2(x1,x2):
    #Calculates the mean differences between two time series x1 and x2. Does it for
    #each row of arrays x1 and x2.
    
    n_rows = np.shape(x1)[0]
    md_out = np.empty(n_rows)
    
    for ii in range(0,n_rows):
        md_out[ii] = md(x1[ii,:],x2[ii,:])
    
    return md_out

def find_optima(x, alpha=6):
    #Finds minima and maxima of a discrete time series dataset. Finds all points
    #that are larger or smaller than all neighbouring values in a window of width
    #alpha.

    L = len(x)
    max_val = np.array([])
    min_val = np.array([])
    max_ii = np.array([])
    min_ii = np.array([])
    
    x_max = np.array(x)
    x_min = np.array(x)
    
    x_max[np.isnan(x)] = -10000
    x_min[np.isnan(x)] = 10000
    
    ii = alpha
    
    while ii in range(alpha,L-alpha-1):
        #Find maxima
        if x[ii] >= np.max(x_max[ii-alpha:ii+alpha]):
            max_val = np.append(max_val,x[ii])
            max_ii = np.append(max_ii,ii)
            ii = ii + alpha - 1
        #Find Minima
        elif x[ii] <= min(x_min[ii-alpha:ii+alpha]):
            min_val = np.append(min_val,x[ii])
            min_ii = np.append(min_ii,ii)
            ii = ii + alpha - 1
        else:
            ii = ii + 1
            
    return max_val, max_ii, min_val, min_ii

def find_optima2(X, alpha=6):
    #2D Find optima. Uses rows as time series.
    
    n_rows, n_cols = np.shape(X)
    n_peaks = n_cols
    max_val = np.empty((n_rows,n_peaks))
    max_ii  = np.empty((n_rows,n_peaks))
    min_val = np.empty((n_rows,n_peaks))
    min_ii  = np.empty((n_rows,n_peaks))
    
    max_val[:] = np.nan
    max_ii[:] = np.nan
    min_val[:] = np.nan
    min_ii[:] = np.nan

    for ii in range(0,n_rows):
        x_tmp = X[ii,:]
        max_val0, max_ii0, min_val0, min_ii0 = find_optima(x_tmp,alpha)
        Lmax = len(max_val0)
        Lmin = len(min_val0)
        max_val[ii,0:Lmax] = max_val0
        max_ii[ii,0:Lmax] = max_ii0
        min_val[ii,0:Lmin] = min_val0
        min_ii[ii,0:Lmin] = min_ii0
        
    return max_val, max_ii, min_val, min_ii

def find_optima_in_window(x, max_ind, min_ind, alpha=3):
    #Finds optima in a range defined by an alpha window either side of indices
    #in ind.
    
    Lmax = len(max_ind)
    Lmin = len(min_ind)
    
    max_val = np.empty(Lmax)
    min_val = np.empty(Lmin)
    max_val[:] = np.nan
    min_val[:] = np.nan
    
    max_ii = np.empty(Lmax)
    min_ii = np.empty(Lmin)
    max_ii[:] = np.nan
    min_ii[:] = np.nan
    
    for ii in range(0,Lmax):
        intind = max_ind[ii] 
        if not np.isnan(intind):
            intind = int(intind)
            max_val[ii] = np.max(x[intind-alpha:intind+alpha])
            max_ii[ii] = np.argmax(x[intind-alpha:intind+alpha]) - alpha
            
    for ii in range(0,Lmin):
        intind = min_ind[ii]
        if not np.isnan(intind):
            intind = int(intind)
            min_val[ii] = np.min(x[intind-alpha:intind+alpha])
            min_ii[ii] = np.argmin(x[intind-alpha:intind+alpha]) - alpha
            
    max_ii = max_ind + max_ii
    min_ii = min_ind + min_ii
    
    return max_val, max_ii, min_val, min_ii

def find_optima_in_window2(X,max_ind,min_ind,alpha=3):
    
    n_rows, n_cols = np.shape(X)
    n_peaks = n_cols
    max_val = np.empty((n_rows,n_peaks))
    max_ii  = np.empty((n_rows,n_peaks))
    min_val = np.empty((n_rows,n_peaks))
    min_ii  = np.empty((n_rows,n_peaks))
    
    max_val[:] = np.nan
    max_ii[:] = np.nan
    min_val[:] = np.nan
    min_ii[:] = np.nan

    for ii in range(0,n_rows):
        max_val0, max_ii0, min_val0, min_ii0 = find_optima_in_window(X[ii,:], max_ind[ii,:], min_ind[ii,:], alpha)
        Lmax = len(max_val0)
        Lmin = len(min_val0)
        max_val[ii,0:Lmax] = max_val0
        min_val[ii,0:Lmin] = min_val0
        max_ii[ii,0:Lmax] = max_ii0
        min_ii[ii,0:Lmin] = min_ii0
        
    return max_val, max_ii, min_val, min_ii
    
    return

def remove_duplicates(a):
    return [ii for ii in a if a.count(ii) == 1]

