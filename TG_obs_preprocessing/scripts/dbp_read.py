"""
This module is for functions that read in specific filetypes.

KEEP TO THE ANGLE CONVENTION -180 -> 180


Changes:
2019-11-21: (dbyrne)    Initial version.
2025-06-04: (jelt)  Ported to NEMO_validation, to be used for model assessment, with minor edits
2025-06-04: (jelt)  Added read_ticon3 function to read TICON-3 data
"""

import numpy as np
import pandas as pd
import dbp_general as dbg
from netCDF4 import Dataset as ds
import matplotlib.pyplot as plt
import dbp_gtm as gtm
from os import listdir
import GTM_assimilation as ass

###
###*****GENERAL READING FUNCTIONS*****
###

def write_nc(filename, varlist, varnames, dimlist, dim_dimensions):
    '''
    # Writes a list of variables to new netcdf file. 
    #
    # filename       :: directory and filename of new netcdf file.
    # varlist        :: list of 
    # dimlist        ::
    # dim_dimensions ::
    '''
    return

def read_nc(filename, varlist, apply_fillvalues=True, custom_fillvalues = [],
            rot_int=False):
    '''
    # Reads an arbitrary number of variables from an NetCDF file and does some 
    # processing. 
    #
 IN # filename           :: netcdf file location.
 IN # varlist            :: list of variables to read from the file.
 IN # apply_fillvalues   :: read and convert fillvalues (from file) to np.nan.
 IN # custom_fillevalues :: convert regions of custom values to np.nan (list).
 IN # rot_int            :: rotate so lons/lats are intuitively rotated. This 
    #                       means that using plt.imshow will show the array
    #                       with North at the top and East on the right.
    #
OUT # var                :: A list of variables read from file, in same order
    #                       as specified in varlist
    '''
    
    n_vars = len(varlist)
    ncid = ds(filename)
    var = [np.array(ncid.variables[ii][:,:]) for ii in varlist]
    
    if apply_fillvalues:
        fills = [np.array(ncid.variables[ii]._FillValue) for ii in varlist]
        for ii in range(0,n_vars):
            var[ii] = dbg.apply_fill(var[ii],fills[ii])
            
    if len(custom_fillvalues) > 0:
        print('Applying Custom Fill Values')
        for ii in range(0,n_vars):
            var[ii] = dbg.apply_fill(var[ii], custom_fillvalues)
        
    if rot_int:
        print('Rotation not yet implemented')
    
    ncid.close()
    return var

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
    
def read_obs_nc(const, fn_dir=None):
    
    if fn_dir is None:
        fn_dir = '/projectsa/NOCglobaltide/data/obs/for_DA/'
    
    fn = fn_dir + 'obs_' + const.upper() + '.nc'
    
    ncid = ds(fn)
    lon = ncid.variables['longitude'][:]
    lat = ncid.variables['latitude'][:]
    z1  = ncid.variables['z1'][:]
    z2  = ncid.variables['z2'][:]
    try:
        obs_id1 = ncid.variables['obs_id1'][:]
        obs_id2 = ncid.variables['obs_id2'][:]
    except:
        obs_id1 = 0
        obs_id2 = 0
    a,g = dbg.cart2polar(z1,z2)
    
    return lon, lat, z1, z2, a, g, obs_id1, obs_id2

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
        for variant in (cc.lower(), cc.upper()):
            try:
                ncid = ds(file_dir + variant + '.nc')
                break
            except FileNotFoundError:
                continue

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
    z1, z2 = dbg.polar2cart(fes_a, fes_g)
    return fes_a, fes_g, fes_lon, fes_lat, fes_mask, z1, z2

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
        nemo_z1 = dbg.apply_mask(nemo_z1,nemo_mask)
        nemo_z2 = dbg.apply_mask(nemo_z2,nemo_mask)
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
        
    nemo_a, nemo_g = dbg.cart2polar(nemo_z1, nemo_z2)
    nemo_g = -nemo_g
    nemo_a = np.squeeze(nemo_a); nemo_g = np.squeeze(nemo_g)
    nemo_z1 = np.squeeze(nemo_z1); nemo_z2 = np.squeeze(nemo_z2)
    
    return nemo_a, nemo_g, nemo_lon, nemo_lat, nemo_mask, nemo_z1, nemo_z2

def read_mdp_tg_harm(file_read, use_datatypes, const):
    '''
    # For reading TG Harmonic data from MDP file. There are three datatypes, 
    # 1,2 and 3. A constituent dictionary must be supplied since the file is 
    # indexed by Doodson index. Code is a bit of a mess and has been
    # added to lazily on a number of occasions.
    '''
    n_mdp_const = 40 # Number of constituents for each location in file
    const_dict = {'SA':1,'SSA':2,'MM':3,'MSF':4,'MF':5,'Q1':8,'O1':10,
                  'P1':15,'S1':16,'K1':17,'J1':21,'2N2':26,'MU2':27,
                  'N2':28,'NU2':29,'M2':31,'MKS2':32,'LA2':33,'L2':34,
                  'T2':35,'S2':36,'R2':37,'K2':38,'M3':43,'MN4':47,'M4':48,
                  'MS4':50,'S4':52,'M6':55,'M8':89}
    n_const = len(const)
    const_id = np.zeros(n_const)*np.nan
    for ii in range(0,n_const):
        const_id[ii] = const_dict[const[ii]]
    const_dict = dict(zip(const, const_id))
    
    #Unpack constituent dictionary and do other dictionary stuff
    out_const = list(const_dict.keys())
    out_const_id = list(const_dict.values())
    n_out_const = len(out_const)
    out_tide_dict = dict( zip( out_const_id, np.arange(0, n_out_const) ) )
    
    #Open and read all of file.
    f_id = open(file_read)
    all_lines = f_id.readlines()
    
    #Determine number of lines/locations
    n_lines = len(all_lines)
    
    #Initialise datasets.
    mdp_data_type  = []
    mdp_lat     = []
    mdp_lon     = []
    mdp_tg_z0      = []
    mdp_port       = []
    indices_tmp    = []
    mdp_amplitude = np.zeros((n_out_const, n_lines))
    mdp_phase     = np.zeros((n_out_const, n_lines))
    mdp_amplitude[:] = np.nan
    mdp_phase[:] = np.nan
    datatype_list = []
    
    #Read  through each line, putting data into the expected arrays.
    line_ii = -1
    for line in all_lines:
        
        line_ii = line_ii + 1
        
        line0                  = line[0:30].split()
        mdp_data_type_tmp      = float(line0[0])
        
        if mdp_data_type_tmp in use_datatypes:
            mdp_data_type.append( mdp_data_type_tmp )
            mdp_lat.append( float(line0[1]) )
            mdp_lon.append( float(line0[2]) )
            mdp_tg_z0.append( float(line0[3]) )
            mdp_port.append( str(line[29:62]).strip() )
            indices_tmp.append( line_ii )
            datatype_list.append( mdp_data_type_tmp )
        
            line0 = line[63:].split()
            for const_ii in range(0,n_mdp_const):
                mdp_tide_id = float(line0[const_ii*3+2])
                
                if mdp_tide_id in out_const_id:
                    out_ii = out_tide_dict[mdp_tide_id]
                    mdp_amplitude[out_ii, line_ii] = float(line0[const_ii*3])
                    mdp_phase[out_ii, line_ii]     = float(line0[const_ii*3+1])
    
    #Convert all to numpy arrays
    mdp_data_type = np.array(mdp_data_type)
    mdp_port      = np.array(mdp_port)
    mdp_lat    = np.array(mdp_lat) 
    mdp_lon    = np.array(mdp_lon)
    mdp_tg_z0     = np.array(mdp_tg_z0)
    
    mdp_a = mdp_amplitude[:,indices_tmp]
    mdp_g = mdp_phase[:,indices_tmp]
    
    #Change interval phase lies in
    mdp_g = mdp_g - 360*(mdp_g>180)
    f_id.close()
    
    mdp_z1, mdp_z2 = dbg.polar2cart(mdp_a, mdp_g)
    
    return mdp_a, mdp_g, mdp_lon, mdp_lat, mdp_z1, mdp_z2

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
        z1, z2 = dbg.polar2cart(amp, pha, degrees=True)
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
    
    z1, z2 = dbg.polar2cart(a,g)
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
    z1, z2 = dbg.polar2cart(a_out, g_out)
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
    file_list = listdir(file_tide_dir)
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
    
    z1, z2 = dbg.polar2cart(a, g)
    
    return a, g, lon, lat, z1, z2

def read_gesla2_file(filename):
    #Reads all information from a standard GESLA2 text file.
    
    f = open(filename)    
    f.readline()
    
    #Read location name
    line = f.readline()
    loc_name = line[12:].strip()
    
    file_skiplines(f,2)
    
    #Read Longitude and Latitude
    line = f.readline().split()
    latitude = line[2]
    line = f.readline().split()
    longitude = line[2]
    
    file_skiplines(f,7)
    
    #Read null_value
    line = f.readline().split()
    
    #skip to data
    file_skiplines(f,18)
    
    date = []
    hour = []
    z    = []
    qc   = []
    
    for line in f.readlines():
        line = line.split()
        date.append(line[0])
        hour.append(line[1])
        z.append(float(line[2]))
        qc.append(float(line[3]))
        
    f.close()
    
    n_obs = len(date)
    time = np.array((1,n_obs))
        
    z = np.array(z)
    qc = np.array(qc)
        
    tmp = np.where(qc!=1)
    z[tmp] = np.nan
    
    for tt in range(0,n_obs):
        yr = date[tt][0:3]
        mm = date[tt][5:6]
        dd = date[tt][8:9]
        hh = hour[tt][11:12]
        mmin = hour[tt][14:15]
        ss = hour[tt][17:18]
        print(yr)
        time[tt] = dbg.date_parser(yr,mm,dd,hh,mmin,ss)
    
    return (loc_name, latitude, longitude, time, z, qc)
    
    
    
