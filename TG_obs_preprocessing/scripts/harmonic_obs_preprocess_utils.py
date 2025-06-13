# harmonic_obs_preprocess_utils.py

"""
Code adapted from: https://github.com/NOC-MSM/global_tide
@author: dbyrne

David Byrne, Jeff Polton & Colin Bell (2021): Creation of a global tide analysis
dataset: Application of NEMO and an offline objective analysis scheme, Journal of Operational
Oceanography, https://doi.org/10.1080/1755876X.2021.2000249


"""

import numpy as np
from netCDF4 import Dataset as ds
import pandas as pd
import os
import time




###
###*****GENERAL READING FUNCTIONS*****
###

def read_obs(cc, dir='/projectsa/NOCglobaltide/data/obs/'):
    '''
    Filthy dirty routine for reading in all of my obs
    '''
    
    fn_psmsl_tide_dir = dir + 'psmsl_bpr/tides/'
    fn_psmsl_locs_dir = dir + 'psmsl_bpr/data/'
    fn_iapso_data     = dir + 'IAPSO/iapso_pelagic_data.txt'
    #fn_mdp_data       = dir + 'TG/Coastal_Points_ALL.txt'
    fn_gloup_head     = dir + 'GLOUP/gloup_headers.dat'
    fn_gloup_harm     = dir + 'GLOUP/gloup_const_all.dat'
    fn_ticon3_data    = dir + 'TICON_3/TICON_3.txt'
    mdp_data_types=[3]
    
    # # # # # READ PSMSL HARMONICS # # # # #
    tmp = read_psmsl_bpr_harm(fn_psmsl_tide_dir, fn_psmsl_locs_dir, cc)
    psmsl_a = tmp[0]; psmsl_g = tmp[1]; psmsl_lon = tmp[2] + 360; 
    psmsl_lat = tmp[3]; psmsl_z1 = tmp[4]; psmsl_z2 = tmp[5]
    psmsl_lon[psmsl_lon>433] = psmsl_lon[psmsl_lon>433] - 360
    if(1):
        # # # # # READ IAPSO HARMONICS # # # # # 
        tmp = read_iapso_bpr_harm(fn_iapso_data, cc)
        iapso_a = tmp[0]; iapso_g = tmp[1]; iapso_lon = tmp[2]+360; 
        iapso_lat = tmp[3]; iapso_z1 = tmp[4]; iapso_z2 = tmp[5]
        iapso_lon[iapso_lon>433] = iapso_lon[iapso_lon>433] - 360
        # # # # #  READ GLOUP HARMONICS # # # # # 
        tmp = read_gloup_bpr_harm(fn_gloup_harm, fn_gloup_head, cc)
        gloup_a = tmp[0]; gloup_g = tmp[1]; gloup_lon = tmp[2]+360; 
        gloup_lat = tmp[3]; gloup_z1 = tmp[4]; gloup_z2 = tmp[5]
        gloup_lon[gloup_lon>433] = gloup_lon[gloup_lon>433] - 360
    # # # # # READ MDP HARMONICS # # # # #
    #tmp = dbr.read_mdp_tg_harm(fn_mdp_data, mdp_data_types, cc)
    #mdp_a = tmp[0]; mdp_g = tmp[1]; mdp_lon = tmp[2]+360; mdp_lat = tmp[3]
    #mdp_z1 = tmp[4]; mdp_z2 = tmp[5]
    #mdp_lon[mdp_lon>433] = mdp_lon[mdp_lon>433] - 360
    # # # # # READ TICON-3 HARMONICS # # # # #
    tmp = read_ticon3(fn_ticon3_data, cc)
    ticon_a = tmp[0]  # amplitude
    ticon_g = tmp[1]  # phase
    ticon_lon = tmp[2]  # longitude
    ticon_lat = tmp[3]  # latitude
    ticon_z1 = tmp[4]  # z1 component
    ticon_z2 = tmp[5]  # z2 component

    lon = np.concatenate((iapso_lon, psmsl_lon, gloup_lon, ticon_lon))
    lat = np.concatenate((iapso_lat, psmsl_lat, gloup_lat, ticon_lat))
    z1 = np.concatenate((iapso_z1, psmsl_z1, gloup_z1, ticon_z1), 1)
    z2 = np.concatenate((iapso_z2, psmsl_z2, gloup_z2, ticon_z2), 1)
    a = np.concatenate((iapso_a, psmsl_a, gloup_a, ticon_a), 1)
    g = np.concatenate((iapso_g, psmsl_g, gloup_g, ticon_g), 1)
    z1 = np.squeeze(z1);z2 = np.squeeze(z2);a=np.squeeze(a);g=np.squeeze(g) 
    
    obs_id = np.zeros(len(lon))
    curind = 0; nexind = curind + len(psmsl_lon)
    obs_id[ :len(psmsl_lon) ] = 0; 
    curind = nexind; nexind = curind + len(iapso_lon)
    obs_id[ curind:nexind ] = 1 
    curind = nexind; nexind = curind + len(gloup_lon)
    obs_id[ curind:nexind ] = 2
    curind = nexind; nexind = curind + len(ticon_lon)
    obs_id[ curind:nexind ] = 3
    
    return lon, lat, z1, z2, a, g, obs_id


def file_skiplines(f,num):
    '''
    # Skips num number of lines in open file f
    '''
    for iii in range(0,num):
        f.readline()
    return

def process_nemo_nc(X,fillvalues=[1e40],flip=True):
    '''
    # Very basic processing of NEMO output .nc data. Flips the data to a 
    # more intuitive orientation and applies any fill values. This is done
    # AFTER reading the data from NC file.
    #
 IN # X                  :: Array to be processed
 IN # fillvalues         :: A list of fill values to be applied.
 IN # flip               :: To flip or not to flip
    #
OUT # var                :: A list of variables read from file, in same order
    #                       as specified in varlist
    '''

    Xout = np.array(X)
    Xout = np.squeeze(Xout)
    
    if flip:
        Xout = np.flip(Xout,0)
    
    for ii in range(0,len(fillvalues)):
        ff = (Xout==fillvalues[ii])
        Xout[ff] = np.nan
    
    return Xout


###
###*****SPECIFIC READING FUNCTIONS*****###
###
    
def read_fes_2D_harm(file_dir, const, istride=4, jstride=8):
    """
    Reads in 2D harmonic files from a directory of FES2014 files.
    Save amplitude in metres, phase, longitude and latitude in degrees.
    """
    #Cast constituents to lower case
    const = [const[ii].lower() for ii in range(0,len(const))]
    
    #Reads in 2D harmonic files from a directory of FES2014 files. 
    fes_a = []
    fes_g = []
    
    #Read in constituents one at a time from FES directory
    #Loop over constituents
    cc_ii = -1
    for cc in const:
        cc_ii = cc_ii + 1
        #Define filename and load
        ncid = ds(file_dir + cc + '.nc')
        #Read in fill value from file
        fill = ncid.variables['amplitude']._FillValue
        #Read and process the relevant slice of nc file
        amp_tmp = np.array(ncid.variables['amplitude'][:][::istride,::jstride])
        pha_tmp = np.array(ncid.variables['phase'][:][::istride,::jstride])
        amp_tmp = process_nemo_nc(amp_tmp, [fill])
        pha_tmp = process_nemo_nc(pha_tmp, [fill])
        #Append constituent to list and convert to metres
        fes_a.append(amp_tmp/100)
        fes_g.append(pha_tmp)
        
    #Define lat/lon arrays
    fes_lat = ncid.variables['lat_bnds'][:][::istride]
    fes_lon = ncid.variables['lon_bnds'][:][::jstride]
    fes_lon2, fes_lat2 = np.meshgrid(fes_lon[:,1],fes_lat[:,0])
    
    fes_lon = np.flip(fes_lon2,0)
    fes_lat = np.flip(fes_lat2,0)
    
    #Estimate FES_mask (1 = land)
    fes_mask = np.array(ncid.variables['amplitude'][:][::istride,::jstride])
    fes_mask = process_nemo_nc(fes_mask)
    fes_mask = fes_mask == fill
    
    fes_a = np.array(fes_a); fes_g = np.array(fes_g); fes_lon = np.array(fes_lon)
    fes_lat = np.array(fes_lat); fes_mask = np.array(fes_mask)
    z1, z2 = polar2cart(fes_a, fes_g)
    return fes_a, fes_g, fes_lon, fes_lat, fes_mask, z1, z2

def read_nemo_2D_grid(file_dom = '', fix_lons=True, 
                      istride = 1, jstride = 1):
    '''
    # Reads and returns lat, lon and mask.
    In progress. Does not replicate the required functionality of read_nemo_2D_harm
    '''

    ncdom = ds(file_dom)
    top_level = ncdom.variables['top_level'][:,:]
    top_level = top_level.astype(float)
    top_level = process_nemo_nc(top_level)
    top_level = top_level[::istride, ::jstride]
    nemo_mask = top_level==0
        
    #Read in latitude, longitude variables from domain file at t-points
    nemo_lat = ncdom.variables['gphit'][:,:]
    nemo_lon = ncdom.variables['glamt'][:,:]

    nemo_lat = nemo_lat[:,::istride, ::jstride]
    nemo_lon = nemo_lon[:,::istride, ::jstride]
        
    #close files
    ncdom.close()
    
    nemo_lon = np.array(nemo_lon)
    nemo_lat = np.array(nemo_lat)

    # Make lons be between -180 and 180 (can't remember where this came from)
    if(0):# fix_lons:
        fixed_lons = nemo_lon.copy()
        for i, start in enumerate(np.argmax(np.abs(np.diff(nemo_lon)) > 180, axis=1)):
            fixed_lons[i, start+1:] += 360        
        nemo_lon = fixed_lons
        
    return nemo_lon, nemo_lat, nemo_mask

def read_nemo_2D_harm(file_harm, const, file_dom = '', fix_lons=True, 
                      istride = 1, jstride = 1, var='z', read_dom = True):
    '''
    # Reads harmonic ssh constituents (z1, z2) from a NEMO output file, 
    # does some processing and converts to amplitude and phase.
    '''
    const_flip = ['K1', 'J1'] # These must be multiplied by -1 to fit obs
    
    cc = 0
    varcopy = var
    if var=='z':
        var = ''
    else:
        var = '_'+var
    
    ncid = ds(file_harm)
    nemo_z1_out = []
    nemo_z2_out = []
    
    if read_dom:
        ncdom = ds(file_dom)
        top_level = ncdom.variables['top_level'][:,:]
        top_level = top_level.astype(float)
        top_level = process_nemo_nc(top_level)
        top_level = top_level[::istride, ::jstride]
        nemo_mask = top_level==0
    
    print('Reading file ' + file_harm)
    #Loop over constituents
    for ii in const:
        #Read in z1, z2 values for current constituent
        nemo_z1 = ncid.variables[const[cc]+'x'+var][:,:][::istride, ::jstride]
        nemo_z2 = ncid.variables[const[cc]+'y'+var][:,:][::istride, ::jstride]
        #Determine the fill value for the constituent and sent to process
        try:
            fill = ncid.variables['M2x']._FillValue
        except:
            fill = ncid.variables['SSAx']._FillValue
        nemo_z1 = process_nemo_nc(nemo_z1,[fill])
        nemo_z2 = process_nemo_nc(nemo_z2,[fill])
        
        #TEMPORARY
        nemo_mask = nemo_z1 == 0
        
        #Apply nan mask based on landmask defined earlier
        nemo_z1 = apply_mask(nemo_z1,nemo_mask)
        nemo_z2 = apply_mask(nemo_z2,nemo_mask)
        # If necessary, flip z1 and z2
        if ii in const_flip:
            nemo_z1 = -nemo_z1
            nemo_z2 = -nemo_z2
        #Convert phases to interval (-180,180)
        #nemo_g_tmp = nemo_g_tmp - 360*(nemo_g_tmp>180)
        nemo_z1_out.append(nemo_z1)
        nemo_z2_out.append(nemo_z2) 
        cc = cc+1
    
    #Read in latitude, longitude variables from domain file
    if read_dom:
        if varcopy == 'z':
            nemo_lat = ncdom.variables['gphit'][:,:]
            nemo_lon = ncdom.variables['glamt'][:,:]
        elif varcopy == 'u':
            nemo_lat = ncdom.variables['gphiu'][:,:]
            nemo_lon = ncdom.variables['glamu'][:,:]
        elif varcopy == 'v':
            nemo_lat = ncdom.variables['gphiv'][:,:]
            nemo_lon = ncdom.variables['glamv'][:,:]
    
        nemo_lat = nemo_lat[:,::istride, ::jstride]
        nemo_lon = nemo_lon[:,::istride, ::jstride]
        #Process latitude longitude variables
        nemo_lat = process_nemo_nc(nemo_lat,[fill])
        nemo_lon = process_nemo_nc(nemo_lon,[fill])
        
    #close files
    ncid.close()
    if read_dom: ncdom.close()
    
    nemo_z1 = np.array(nemo_z1_out); nemo_z2 = np.array(nemo_z2_out)
    nemo_mask = np.array(nemo_mask)
    
    if read_dom:
        nemo_lon = np.array(nemo_lon); nemo_lat = np.array(nemo_lat); 
    else:
        nemo_lon = np.zeros(nemo_z1_out[0].shape)
        nemo_lat = np.zeros(nemo_z1_out[0].shape)
    
    # Make lons be between -180 and 180 (can't remember where this came from)
    if fix_lons:
        fixed_lons = nemo_lon.copy()
        for i, start in enumerate(np.argmax(np.abs(np.diff(nemo_lon)) > 180, axis=1)):
            fixed_lons[i, start+1:] += 360        
        nemo_lon = fixed_lons
        
    nemo_a, nemo_g = cart2polar(nemo_z1, nemo_z2)
    nemo_g = -nemo_g
    nemo_a = np.squeeze(nemo_a); nemo_g = np.squeeze(nemo_g)
    nemo_z1 = np.squeeze(nemo_z1); nemo_z2 = np.squeeze(nemo_z2)
    
    return nemo_a, nemo_g, nemo_lon, nemo_lat, nemo_mask, nemo_z1, nemo_z2

def read_ticon3(fn_ticon3_data, cc, gauge_type='Coastal'):
    """
    Reads the TICON-3 data file and returns the amplitude, phase, longitude, latitude,
     z1 and z2 components for a specified constituent and gauge type.

    If amplitude/phase are flagged or not present, they will be set to NaN.
    Latitude and longitude will be set to zero for this request.

    Save amplitude in metres, phase, longitude and latitude in degrees.

    Data source:
    Hart-Davis, Michael G; Dettmering, Denise; Seitz, Florian (2022): TICON-3: Tidal Constants based on GESLA-3 sea-level records from globally distributed tide gauges including gauge type information (data) [dataset]. PANGAEA, https://doi.org/10.1594/PANGAEA.951610

    Changes:
    2025-06-04: Created (jelt)
    """
    

    # load data
    #fn_ticon3_data = "/Users/jelt/Downloads/TICON_3/TICON_3.txt" # TICON-3 file
    df = pd.read_csv(fn_ticon3_data, sep='\t', header=None)
    
    # initialize variables
    lat = 0; lon = 0; amp = np.nan; pha = np.nan; z1 = np.nan; z2 = np.nan

    try:
        # select a constituent
        if isinstance(cc, list) and len(cc) == 1:
            cc = cc[0]
        indx_cons = np.where(df[2] == cc)
        df = df.iloc[indx_cons]
        # select tide gauge type [options are: River, Lake and Coastal]
        indx_type = np.where(df[13] == gauge_type)
        df = df.iloc[indx_type]
        lon = np.array(df[1])  # assign longitude
        lat = np.array(df[0]) # assign latitude
        amp = np.array(df[3])/100.; amp = amp.reshape(1, len(amp)) # assign amplitude
        pha = np.array(df[4]); pha = pha.reshape(1, len(pha)) # assign phase
        #con = np.array(df[2]); con = con.reshape(1, len(con) # assign consituent name
        # Also change interval phase lies in from 0 -> 360 to -180 -> 180.
        pha[pha>180] = pha[pha>180] - 360
        # Convert amplitudes and phases to z1 and z2
        z1, z2 = polar2cart(amp, pha, degrees=True)
    except Exception as e:
        pass
    return amp, pha, lon, lat, z1, z2

def read_gloup_bpr_harm(fn_harm, fn_head, const):
    '''
    # Reads in data from GLOUP BPR textfile.
    #
    # If amplitude/phase are flagged or not present, they will be set to NaN.
    # Latitude and longitude will still be saved for this point however.
    #
    # Saves amplitude in metres and phase, longitude and latitude in degrees.

    # Author: dbyrne | Version 1.0 (21/11/2019)
    '''
    
    f_harm = open(fn_harm)
    f_head = open(fn_head)
    n_locs = 243
    n_const = len(const)
    const_dict = dict(zip(const,np.arange(0,n_const)))
    
    a = np.zeros((n_const, n_locs))*np.nan
    g = np.zeros((n_const, n_locs))*np.nan
    lat = np.zeros(n_locs)*np.nan
    lon = np.zeros(n_locs)*np.nan
    n_ana_const = np.zeros(n_locs)*np.nan
    
    file_skiplines(f_head,11)
    keep_index = []
    
    for ii in range(0,n_locs):
        lat[ii] = float( f_head.readline()[:9] )
        lon[ii] = float( f_head.readline()[:9] )
        file_skiplines(f_head,25)
        
        keep_index.append(ii)
        loop = True
        while loop:
            line = f_harm.readline()[1:6]
            if line == 'MAJOR': 
                loop = False
                file_skiplines(f_harm, 4)
        
        loop = True
        loop_count = 0
        while loop:
            line = f_harm.readline()
            if line[1:6] == 'INPUT' or line[:3] == 'END':
                loop = False
            else:
                const_tmp = line[5:10].strip()
                loop_count = loop_count + 1
                if const_tmp in const:
                    #print(const_tmp)
                    index = const_dict[const_tmp]
                    a[index, ii] = float( line[27:39] )/100
                    g[index, ii] = float( line[39:50] )
        n_ana_const[ii] = loop_count
                
    f_harm.close()
    f_head.close()
    a = a[:,keep_index]; g = g[:,keep_index];
    lon = lon[keep_index]; lat = lat[keep_index];
    
    z1, z2 = polar2cart(a,g)
    g[g>180] = g[g>180] - 360
    
    return a, g, lon, lat, z1, z2, n_ana_const

def read_iapso_bpr_harm(file_read, const, include_flagged=False):
    '''
    # Reads in data from IAPSO BPR textfile. Each location has the same 
    # constituent variable order, given in gloup_dict. The input const is the 
    # desired constituent order for the output. If constituent is not contained
    # inside the GLOUP file then that part of the array will be a nan column.
    # Include_flagged currently has no real functionality but should be kept as
    # True.
    #
    # If amplitude/phase are flagged or not present, they will be set to NaN.
    # Latitude and longitude will still be saved for this point however.
    #
    # Author: dbyrne | Version 1.0 (21/11/2019)
    '''
    
    # Locations of each constituent in GLOUP file.
    gloup_dict = {'Q1':0,'O1':1,'P1':2,'K1':3,'N2':4,'M2':5,'S2':6,'K2':7}
    n_gloup = len(gloup_dict)
    
    #Read in all lines, determine number of locations and initialise arrays
    f = open(file_read)
    all_lines = f.readlines()
    n_loc = int(len(all_lines)/3)
    a = np.zeros((n_gloup,n_loc)); g = np.zeros((n_gloup,n_loc))
    lat = np.zeros(n_loc); lon = np.zeros(n_loc)
    
    for ii in range(0,n_loc): #Loop over locations
        
        tmpline = all_lines[ii*3].split()
        if tmpline[1] != '*' or include_flagged:
            
            # Phantom star dummy index for modifying index if data flag present
            phs = 0
            if tmpline[1] == '*':
                phs = 1 
            
            # Read in latitude/longitudes for location and *-1 if necessary
            lat[ii] = float(tmpline[1+phs]) + float(tmpline[2+phs])/60
            if tmpline[3+phs].strip() == 'S':
                lat[ii] = -1*lat[ii]
            lon[ii] = float(tmpline[4+phs]) + float(tmpline[5+phs])/60
            if tmpline[6+phs].strip() == 'W':
                lon[ii] = -1*lon[ii]
            
            # Read in amplitude/phase from first data line
            # Also change interval phase lies in from 0 -> 360 to -180 -> 180.
            a[:,ii] = [float(ii) for ii in all_lines[ii*3+1][1:].split()]
            # Next line
            g[:,ii] = [float(ii) for ii in all_lines[ii*3+2][1:].split()]
        else:
            a[:,ii] = np.nan; g[:,ii] = np.nan
    # Missing Values
    a[a==-1] = np.nan
    g[g==-1] = np.nan
    a = a/100
    
    # Map GLOUP data to output arrays
    n_const = len(const)
    a_out = np.zeros((n_const, n_loc))*np.nan
    g_out = np.zeros((n_const, n_loc))*np.nan
    for ii in range(0,n_const):
        cc = const[ii]
        if cc in gloup_dict:
            a_out[ii] = a[gloup_dict[cc]]
            g_out[ii] = g[gloup_dict[cc]]
    
    # Convert amplitudes and phases to z1 and z2
    z1, z2 = polar2cart(a_out, g_out)
    g_out[g_out>180] = g_out[g_out>180] - 360

    return a_out, g_out, lon, lat, z1, z2

def read_psmsl_bpr_harm(file_tide_dir, file_locs_dir, const):
    '''
    # Reads harmonic data from PSMSL tide text files. Will read from all files
    # in the given directory so make sure there's nothing unexpected in there.
    # Fetches location data from corresponding files in the locs_directory.
    #
    # Author: dbyrne | Version 1.0 (21/11/2019)
    '''
    
    n_const = len(const)
    const_dict = dict(zip(const,np.arange(0,n_const)))
    file_list = os.listdir(file_tide_dir)
    n_files = len(file_list)
    
    lat = np.zeros(n_files); lon = np.zeros(n_files)
    a = np.zeros((n_const, n_files))*np.nan 
    g = np.zeros((n_const, n_files))*np.nan
    
    file_ii = 0

    for file_tide in file_list:
        file_pre = file_tide[:-8]
        file_locs = file_pre + '_hrp.txt'
        fid_tide = open(file_tide_dir + file_tide)
        fid_locs = open(file_locs_dir + file_locs)
    
        file_skiplines(fid_locs,4)
        lat[file_ii] = fid_locs.readline().split()[-1]
        lon[file_ii] = fid_locs.readline().split()[-1]
        
        all_lines = fid_tide.readlines()[8:]
        n_lines = len(all_lines)
        all_lines = [all_lines[ii].split() for ii in range(0,n_lines)]
        all_lines = np.array(all_lines)
            
        for ii in range(0,len(all_lines[:,0])):
            const_tmp = all_lines[ii,0]
            if const_tmp in const:
                const_id = const_dict[const_tmp]
                a[const_id,file_ii] = float(all_lines[ii,3])
                g_tmp = float(all_lines[ii,5])
                if g_tmp > 180:
                    g_tmp = g_tmp - 360
                g[const_id,file_ii] = g_tmp 
        
        file_ii = file_ii+1
    
    a = a/100
    
    z1, z2 = polar2cart(a, g)
    
    return a, g, lon, lat, z1, z2


###
###*****WRITING METHODS*****###
###

def obs_harmonic_file_create(fn, n_obs, c):
    '''
    # Creates a new obs harmonic file, overwriting an old one if it exists. The base 
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

def obs_harmonic_file_append(fn, var, varname, varunits = ''):
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


###
###*****WORKFLOW FUNCTIONS*****###
###

def process_obs(y_lon, y_lat, y, obs_id, nemo_lon, nemo_lat, nemo_mask,
                fes_lon, fes_lat, fes_z, fes_mask, grid_obs_rad = 30, thin_obs_rad = 50,
                y_ii = None, y_jj = None, fes_ii = None, fes_jj = None,):
    '''
    Processes obs for assimilation. Multiple stages:
        1. Take only obs who corresponding model points are near.
        2. Compare to FES and take only those obs that are close to FES.
        3. Of remaining obs, take only those from a long enough analysis.
        4. Then run obs through superobs routine to thin obs.
    For the GTM assimilation scheme, this is done as a pre processing step
    in the scr_process_obs script, which also saves to .nc file.
    '''
    
    print('Initial #obs: ' +  str(len(y_lon)))
    
    
    # Remove points that are distant from NEMO grid
    if y_ii is None:
        y_ii, y_jj = geo_nn_indices(nemo_lon, nemo_lat, y_lon, y_lat, 
                                        mask = nemo_mask)
    nlon = nemo_lon[y_ii, y_jj]; nlat = nemo_lat[y_ii, y_jj]
    dist = dist_haversine(nlon, nlat, y_lon, y_lat)
    nanii = dist>grid_obs_rad
    y[nanii] = np.nan
    print('Points with bad model grid correspondence: ' + str(sum(nanii)))
    
    # Compare to FES data and take only obs that are within 1/2 error st. dev
    if fes_ii is None:
        fes_ii, fes_jj = geo_nn_indices(fes_lon, fes_lat, y_lon, y_lat, 
                                    mask = fes_mask)
    y_fes = fes_z[fes_ii, fes_jj]
    ef = y - y_fes
    estd = np.nanstd(ef); emean = np.nanmean(ef)
    aef = np.abs(ef - emean)
    nanii = aef > 0.5*estd
    y[nanii] = np.nan
    print('Comparison to FES data made: ' + str( sum(nanii )) + ' removed.')
    
    # Superobs
    y_lon, y_lat, y, obs_id = superobs(y_lon, y_lat, y, crit_dist= thin_obs_rad, 
                               method='delete', obs_id = obs_id)
    print('Superobs processed.')
    
    return y_lon, y_lat, y, obs_id


def superobs(lon, lat, obs, crit_dist = 50, method = 'interp', obs_id=None):
    """
    Thin the observations by comparing distances between all pairs and sequentially removing points 
    until all distances are above crit_dist. The point, from the pair, that is removed is the one closest to any other point.
    
    If method is 'interp', the points (lat,lon,obs) are averaged without recalculating distances
    If 'delete' the point is simply removed.
    
    2021 dbyrne
    13/06/2025 jelt - updated for speed: remove the need to recalculate distances each loop. Actively select point for removal based on distance to other points.
    """
    
    if obs_id is None:
        obs_id = np.zeros(lon.shape)
    
    D = dist_matrix_haversine(lon,lat,lon,lat)
    np.fill_diagonal(D,10000)
    crit = np.nanmin(D)
    
    lon_sup = np.array(lon)
    lat_sup = np.array(lat)
    obs_sup = np.array(obs)
    
    while crit < crit_dist:
        minii = np.unravel_index(np.argmin(D), D.shape)
        ii = minii[0]
        jj = minii[1]
        if np.min(np.delete(D[ii, :], jj)) <= np.min(np.delete(D[:, jj], ii)): # delete the point that is closest to any other
            kk = ii; notkk = jj
        else:
            kk = jj; notkk = ii

        if method == 'interp':
            print('method interp needs to be updated')
            lon_sup[notkk] = np.nanmean([lon_sup[kk], lon_sup[notkk]])
            lat_sup[notkk] = np.nanmean([lat_sup[kk], lat_sup[notkk]])
            obs_sup[notkk] = np.nanmean([obs_sup[kk], obs_sup[notkk]])
            lon_sup = np.delete(lon_sup, kk)
            lat_sup = np.delete(lat_sup, kk)
            obs_sup = np.delete(obs_sup, kk)
        elif method == 'delete':
            lon_sup = np.delete(lon_sup, kk)
            lat_sup = np.delete(lat_sup, kk)
            obs_sup = np.delete(obs_sup, kk)
            obs_id = np.delete(obs_id, kk)
        
        # Instead of recalculating D, just drop row and column kk:
        D = np.delete(D, kk, axis=0)
        D = np.delete(D, kk, axis=1)     
        np.fill_diagonal(D,10000)
        crit = np.min(D)
    
    return lon_sup, lat_sup, obs_sup, obs_id


###
###*****GENERAL USE FUNCTIONS*****###
###


def cart2polar(x, y, degrees = True):
    '''
    # Conversion of cartesian to polar coordinate system
    # Output theta is in radians
    @author: dbyrne
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
    @author: dbyrne
    '''
    if degrees:
        theta = np.deg2rad(theta)
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return(x, y)


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
    @author: dbyrne
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
    @author: dbyrne
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
    @author: dbyrne
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


def reduce_array_using_mask(A, mask):
    '''
    # Removes points from an array A based on values in mask. Input arrays A
    # and mask are 2-dimensional and output is a 1-dimensional reduced array.
    # Reduction is done vertically (along columns first. Hence 'F'). It is the True values
    # in the mask that are removed from A.
    @author: dbyrne
    '''
    A = np.array(A)
    A = A.flatten('F')
    mask = np.array(mask)
    mask = mask.flatten('F')
    A_red = A[ np.where(mask==0) ]
    return A_red


def apply_mask(A,mask):
    # Applies 2D mask to a 2D array A. Replaces locations where mask = 1 with 
    # nan 
    A = np.array(A)
    A[mask] = np.nan
    return A
