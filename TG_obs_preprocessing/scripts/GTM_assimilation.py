#Package of GTM assimilation related routines

import numpy as np
import matplotlib.pyplot as plt
import dbp_general as dbg
import dbp_gtm as gtm
import dbp_plot as dbpp
import numpy.random as rand
import multiprocessing as mp
import dbp_read as dbr
import pandas as pd

##############################
 ####SPECIFIC ROUTINES#######
##############################
 
def compare2FES(fes_z1, fes_z2, fes_lon, fes_lat, y_z1, y_z2, y_lon, y_lat):
    
    y_ii, y_jj = dbg.geo_nn_indices(fes_lon, fes_lat, y_lon, y_lat, 
                                    mask = fes_mask)
    
    return keep_id
 
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
    tmp = dbr.read_psmsl_bpr_harm(fn_psmsl_tide_dir, fn_psmsl_locs_dir, cc)
    psmsl_a = tmp[0]; psmsl_g = tmp[1]; psmsl_lon = tmp[2] + 360; 
    psmsl_lat = tmp[3]; psmsl_z1 = tmp[4]; psmsl_z2 = tmp[5]
    psmsl_lon[psmsl_lon>433] = psmsl_lon[psmsl_lon>433] - 360
    if(1):
        # # # # # READ IAPSO HARMONICS # # # # # 
        tmp = dbr.read_iapso_bpr_harm(fn_iapso_data, cc)
        iapso_a = tmp[0]; iapso_g = tmp[1]; iapso_lon = tmp[2]+360; 
        iapso_lat = tmp[3]; iapso_z1 = tmp[4]; iapso_z2 = tmp[5]
        iapso_lon[iapso_lon>433] = iapso_lon[iapso_lon>433] - 360
        # # # # #  READ GLOUP HARMONICS # # # # # 
        tmp = dbr.read_gloup_bpr_harm(fn_gloup_harm, fn_gloup_head, cc)
        gloup_a = tmp[0]; gloup_g = tmp[1]; gloup_lon = tmp[2]+360; 
        gloup_lat = tmp[3]; gloup_z1 = tmp[4]; gloup_z2 = tmp[5]
        gloup_lon[gloup_lon>433] = gloup_lon[gloup_lon>433] - 360
    # # # # # READ MDP HARMONICS # # # # #
    #tmp = dbr.read_mdp_tg_harm(fn_mdp_data, mdp_data_types, cc)
    #mdp_a = tmp[0]; mdp_g = tmp[1]; mdp_lon = tmp[2]+360; mdp_lat = tmp[3]
    #mdp_z1 = tmp[4]; mdp_z2 = tmp[5]
    #mdp_lon[mdp_lon>433] = mdp_lon[mdp_lon>433] - 360
    # # # # # READ TICON-3 HARMONICS # # # # #
    tmp = dbr.read_ticon3(fn_ticon3_data, cc)
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
 
def read_ensemble_array(dn, fn_list, const, nemo_stride = 1, var='z'):
    '''
    # Reads in an ensemble of harmonic datasets for use in covariance
    # estimation for data assimilation. 
    #
    # dn - directory containing ensemble files
    # fn_list - list of ensemble file names
    # const - constituents to read
    #
    # Author: David Byrne | Version: 1.1 (24/02/2020)
    '''
    
    n_ens = len(fn_list)     # Number of ensemble datasets
    print(str(n_ens) + ' ensemble files..')
    ens_z1 = []
    ens_z2 = []
    # Read model ensemble data, looping through each file
    for ensii in range(0,n_ens):
        fn = dn + fn_list[ensii] + '.nc'
        if var == 'z':
            tmp = dbr.read_nemo_2D_harm(fn, const, istride = nemo_stride, 
                                        jstride= nemo_stride, read_dom=False)
        elif var == 'u':
            tmp = dbr.read_nemo_2D_harm(fn, const, istride = nemo_stride, 
                                        jstride= nemo_stride, read_dom=False,
                                        var = 'u')
        elif var == 'v':
            tmp = dbr.read_nemo_2D_harm(fn, const, istride = nemo_stride, 
                                        jstride= nemo_stride, read_dom=False,
                                        var = 'v')
        tmpz1 = tmp[5]
        tmpz2 = tmp[6]
        tmpz1[np.isnan(tmpz1)] = 0
        tmpz1[np.isnan(tmpz2)] = 0
        ens_z1.append(tmpz1)
        ens_z2.append(tmpz2)
    
    ens_z1 = np.array(ens_z1); ens_z2 = np.array(ens_z2)
    #ens_z1 = ens_z1 - np.nanmean(ens_z1, axis=0)
    #ens_z2 = ens_z2 - np.nanmean(ens_z2, axis=0)
        
    return ens_z1, ens_z2



def GTM_OI_ensemble_obs(obsL, n_obs_ens):
    '''
    Generates indices to remove from an observation dataset. To be used in
    ensemble assimilation runs
    '''
    
    obsL = int(obsL)
    n_obs_ens = int(n_obs_ens)
    X = np.arange(0,obsL)
    rand.shuffle(X)
    remove_obs = np.sort(X[0:n_obs_ens])
    
    return remove_obs

def GTM_OI_by_mask(mod_lon, mod_lat, mod_z, 
                   obs_lon, obs_lat, obs_z, 
                   dom_mask, land_mask, Rvar,
                   Brad = 1000, cov_ens=None,
                   obs_ii = None, obs_jj = None):
    '''
    # Do an EnOI/EnKF assimilation for the Global Tide Model. This is a highly 
    # specific routine and not for general use. B is calculated using SSH 
    # covariance and an arbitrary assimilation decay radius. The model 
    # amplitudes and phases are required for reconstruction of SSH fields. 
    # All distances in km. B is not constructed directly at any point, but 
    # instead BH^T and HBH^T are constructed in order to conserve memory and 
    # allow for assimilation of obs outside of the assimilation domain. cov_alg
    # determines which algorithm to use for calculation of covariance matrix B.
    # If cov_alg == 'ens', then cov_ens must be supplied. This takes the form
    # nz*nr*nc, where nr*nc is the model domain size and nz is the number of 
    # ensemble members.
    #
    # Author: David Byrne | Version 1.1
    '''
    
    # Determine whether obs_ii and obs_jj are known
    obs_loc_known = True
    if obs_ii is None:
        obs_loc_known = False 
    
    # Make copies of input arrays (unsure if needed but J.I.C)
    mod_lon = np.array(mod_lon)
    mod_lat = np.array(mod_lat)
    mod_z = np.array(mod_z)
    obs_lon = np.array(obs_lon)
    obs_lat = np.array(obs_lat)
    land_mask = np.array(land_mask)
    
    # Remove any nans from model data (there shouldn't be any but J.I.C)
    mod_z[np.isnan(mod_z)] = 0
    
    # Extract model variables within assimilation domain defined by dom_mask
    rmin, rmax, cmin, cmax = extract_dom_by_mask(dom_mask)
    mod_lon_dom     = mod_lon[rmin:rmax, cmin:cmax]
    mod_lat_dom     = mod_lat[rmin:rmax, cmin:cmax]
    mod_z_dom       = mod_z[rmin:rmax, cmin:cmax]
    land_mask_dom   = land_mask[rmin:rmax, cmin:cmax]
    
    # Flatten model data (will have 1 row and multiple columns)
    mod_z_domF = mod_z_dom.flatten()
    mod_lon_domF = mod_lon_dom.flatten()
    mod_lat_domF = mod_lat_dom.flatten()
    land_mask_domF = land_mask_dom.flatten()
    
    # Determine various sizes of datasets
    obsL = len(obs_z)
    modL_dom = len(mod_z_domF)
    nr, nc = mod_z.shape
    nr_dom, nc_dom = mod_z_dom.shape
    modL_nonan = len(np.squeeze(np.where(land_mask_domF==0)))
    
    # Only do assimilation heavy lifting if there is an adequate number of obs
    # and model points -- arbitrary
    if obsL >1 and modL_nonan>1:
        print('#Obs: ' + str(obsL) + '| #Mod: ' + str(modL_nonan) )
        # Determine 2D indices of all obs in the global domain (not)
        # assimilation domain)
        if obs_loc_known == False:
            obs_ii, obs_jj = dbg.geo_nn_indices(mod_lon, mod_lat, 
                                                obs_lon, obs_lat,
                                                land_mask)
        
        #Calculate ensemble covariance matrices
        ens_ts_modF = cov_ens[:,rmin:rmax, cmin:cmax].reshape(-1,modL_dom)
        ens_ts_obs  = np.squeeze(cov_ens[:,obs_ii, obs_jj])
        ens_ts_cat = np.concatenate((ens_ts_modF, ens_ts_obs), axis = 1)
        ens_ts_cat = np.transpose(ens_ts_cat)
        C = np.cov(ens_ts_cat)
        C1 = C[modL_dom:, :modL_dom]
        C2 = C[modL_dom:, modL_dom:]
        
        # Get distance matrices for assimilation radius. D1 is the distance
        # between model points from observation locations and D2 is distance 
        # between observation locations.
        mod_obs_lon = mod_lon[obs_ii, obs_jj]
        mod_obs_lat = mod_lat[obs_ii, obs_jj]
        D1 = dbg.dist_matrix_haversine(mod_lon_domF, mod_lat_domF,
                                       mod_obs_lon, mod_obs_lat)
        D2 = dbg.dist_matrix_haversine(mod_obs_lon, mod_obs_lat,
                                       mod_obs_lon, mod_obs_lat)
        
        # Get matrices of gaussian values based on distance matrices
        G1 = dbg.gauss(D1, Brad)
        G2 = dbg.gauss(D2, Brad)
        
        # Cosntruct implicit model error covariance B matrices
        BH = np.transpose(C1*G1)
        HBH = C2*G2 
        
        # Construct a basic diagonal observation covariance matrix, R
        R = np.identity(obsL)*Rvar
        
        # Determine innovations
        innov = obs_z - mod_z[obs_ii, obs_jj]
        
        # Do assimilation
        xa, inc = OI2(mod_z_domF, innov, BH, HBH, R)
            
        # Reinstate land values as nans, reshape output to 2D and place back
        # into global domain
        xa[land_mask_dom.flatten()] = np.nan
        inc[land_mask_dom.flatten()] = np.nan
        xa = xa.reshape((nr_dom, nc_dom))
        inc = inc.reshape((nr_dom, nc_dom))
        #xa = replace_model_var_by_mask(xa, dom_mask)
        #inc = replace_model_var_by_mask(inc, dom_mask)
        
    else:
        # If assimilation did not take place, put nans back into global domain.
        #xa = np.zeros((nr,nc))*np.nan
        #inc = np.zeros((nr,nc))*np.nan
        xa = np.zeros((nr_dom, nc_dom))*np.nan
        inc = np.zeros((nr_dom, nc_dom))*np.nan

    return xa, inc, dom_mask

def GTM_OI_by_index(mod_lon_dom, mod_lat_dom, mod_z_dom, land_mask_dom,
                    mod_lon_obs, mod_lat_obs, mod_z_obs, obs_z, 
                    Rvar, Brad, cov_ens_dom, cov_ens_obs, ii):
    '''
    # Do an EnOI/EnKF assimilation for the Global Tide Model. This is a highly 
    # specific routine and not for general use. B is calculated an ensemble
    # and an arbitrary assimilation decay radius. The model 
    # amplitudes and phases are required for reconstruction of SSH fields. 
    # All distances in km. B is not constructed directly at any point, but 
    # instead BH^T and HBH^T are constructed in order to conserve memory and 
    # allow for assimilation of obs outside of the assimilation domain.
    # Covariance ensemble cov_ens takes the form nz*nr*nc, where nr*nc is the 
    # model domain size and nz is the number of ensemble members.
    #
    # ind_tup = [rmin, cmin, rmax, cmax]
    #
    # Author: David Byrne | Version 1.01 (28/02/2020)
    '''
     
    # Size of all data
    nr_dom, nc_dom = mod_z_dom.shape
    
    # Remove any nans from model data (there shouldn't be any but J.I.C)
    mod_z_dom[np.isnan(mod_z_dom)] = 0
     
    # Flatten model data (will have 1 row and multiple columns)
    mod_z_domF = mod_z_dom.flatten()
    mod_lon_domF = mod_lon_dom.flatten()
    mod_lat_domF = mod_lat_dom.flatten()
    land_mask_domF = land_mask_dom.flatten()
    
    # Determine various sizes of remaining datasets
    obsL = len(obs_z)
    modL_dom = len(mod_z_domF)
    modL_nonan = len(np.squeeze(np.where(land_mask_domF==0)))
    
    # Only do assimilation heavy lifting if there is an adequate number of obs
    # and model points -- arbitrary
    if obsL >1 and modL_nonan>=1:
        try:
                
            #Calculate ensemble covariance matrices
            ens_ts_modF = cov_ens_dom.reshape(-1,modL_dom)
            ens_ts_obs  = cov_ens_obs
            ens_ts_cat = np.concatenate((ens_ts_modF, ens_ts_obs), axis = 1)
            ens_ts_cat = np.transpose(ens_ts_cat)
            C = np.cov(ens_ts_cat)
            C1 = C[modL_dom:, :modL_dom]
            C2 = C[modL_dom:, modL_dom:]
            
            # Get distance matrices for assimilation radius. D1 is the distance
            # between model points from observation locations and D2 is distance 
            # between observation locations.
            D1 = dbg.dist_matrix_haversine(mod_lon_domF, mod_lat_domF,
                                           mod_lon_obs, mod_lat_obs)
            D2 = dbg.dist_matrix_haversine(mod_lon_obs, mod_lat_obs,
                                           mod_lon_obs, mod_lat_obs)
            
            # Get matrices of gaussian values based on distance matrices
            G1 = dbg.gauss(D1, Brad)
            G2 = dbg.gauss(D2, Brad)
            
            # Cosntruct implicit model error covariance B matrices
            BH = np.transpose(C1*G1)
            HBH = C2*G2 
            
            # Construct a basic diagonal observation covariance matrix, R
            R = np.diag(Rvar)
            
            # Determine innovations
            innov = obs_z - mod_z_obs
            
            # Do assimilation
            xa, inc = OI2(mod_z_domF, innov, BH, HBH, R)
                
            # Reinstate land values as nans, reshape output to 2D and place back
            # into global domain
            xa[land_mask_dom.flatten()] = np.nan
            inc[land_mask_dom.flatten()] = np.nan
            xa = xa.reshape((nr_dom, nc_dom))
            inc = inc.reshape((nr_dom, nc_dom))
        except:
            print('Problem with: ' + str(ii))
            xa = np.zeros((nr_dom, nc_dom))*np.nan
            inc = np.zeros((nr_dom, nc_dom))*np.nan
        
    else:
        # If assimilation did not take place, put nans back into global domain.
        xa = np.zeros((nr_dom, nc_dom))*np.nan
        inc = np.zeros((nr_dom, nc_dom))*np.nan

    outstr = str(ii) + ' ::: '
    outstr = outstr + '#Obs: ' + str(obsL) + '| #Mod: ' + str(modL_nonan) 
    print(outstr)

    return xa, inc, ii

def GTM_OI_global(mod_lon, mod_lat, mod_z, land_mask,
                  obs_lon, obs_lat, obs_z, Brad, Rvar, boxsize, 
                  cov_ens, obs_ii = None, obs_jj = None):
    '''
    # For handling the box-wise assimilation into the global GTM domain.
    # This is the original SERIAL version of the routine.
    #
    # Author: David Byrne | Version 1.2
    '''
    # Determine the number of assimilation boxes
    nr, nc = mod_lon.shape
    
    # Initialise arrays for assimilation box storage
    xa_global = np.zeros((nr , nc)) * np.nan
    inc_global = np.zeros((nr , nc)) * np.nan
        
    colrange = np.arange(0, nc, boxsize[1])
    rowrange = np.arange(0, nr, boxsize[0])
    
    # Define list of assimilation domains and number.
    dom_mask = [domain_mask_by_custom(nr, nc, 
                [rowbox, rowbox + boxsize[0] + 1],
                [colbox, colbox + boxsize[1] + 1]) 
                for colbox in colrange for rowbox in rowrange]
    
    dL = len(dom_mask)
    print('Number of domains: ' + str(dL))
    
    # Do the assimilation
    out = [GTM_OI_by_mask(mod_lon, mod_lat, mod_z, obs_lon, obs_lat, obs_z, 
                              dom_maskii, land_mask, Rvar=Rvar, Brad=Brad, 
                              cov_ens=cov_ens, obs_ii = obs_ii, obs_jj = obs_jj)
                              for dom_maskii in dom_mask]
    out = np.array(out)
    xa_global = np.nanmean(out[:,0], axis=0)
    inc_global = np.nanmean(out[:,0], axis=0)
    
    return xa_global, inc_global

def GTM_OI_global_par(mod_lon, mod_lat, mod_z, land_mask,
                  obs_lon, obs_lat, obs_z, Brad, Rvar, boxsize, 
                  cov_ens, obs_ii = None, obs_jj = None, n_proc=24,
                  rowbnds = 'auto'):
    '''
    # For handling the box-wise assimilation into the global GTM domain. Does
    # this using PARALLEL setup. Will do as many assimilation boxes at one
    # time as allowed by the number of processes (n_proc).
    #
    # Author: David Byrne | Version 1.0
    '''
    # Determine the indices and number of assimilation sub-domains.
    # These are stored in dom_index_list = [rmin, rmax, cmin, cmax]
    nr, nc = mod_lon.shape
    if rowbnds == 'auto':
        rowbnds = [0, nr]
    rowrange = np.arange(rowbnds[0], rowbnds[1], boxsize[0])
    colrange = np.arange(0, nc, boxsize[1])
    dom_index_list = np.array( [ [rr, rr + boxsize[0], 
                                  cc, cc + boxsize[1]] 
                                 for rr in rowrange for cc in colrange] )
    dL = dom_index_list.shape[0] 
    for ii in range(0,dL):
        if dom_index_list[ii,1] >= nr: dom_index_list[ii,1] = nr-1
        if dom_index_list[ii,3] >= nc: dom_index_list[ii,3] = nc-1
    print('Domains ready. Number of domains: ' + str(dL) ) 
        
    cov_ens_dom = []
    mod_lon_dom = []
    mod_lat_dom = []
    mod_z_dom = []
    land_mask_dom = []    

    for ii in range(0,dL):
        ind_tmp = dom_index_list[ii]
        rmin = ind_tmp[0]; rmax = ind_tmp[1]; cmin = ind_tmp[2]; cmax = ind_tmp[3]
        cov_ens_dom.append( np.array( cov_ens[:, rmin:rmax, cmin:cmax] ) )
        mod_lon_dom.append( np.array( mod_lon[rmin:rmax, cmin:cmax] ) )
        mod_lat_dom.append( np.array( mod_lat[rmin:rmax, cmin:cmax] ) )
        mod_z_dom.append( np.array( mod_z[rmin:rmax, cmin:cmax] ) )
        land_mask_dom.append( np.array( land_mask[rmin:rmax, cmin:cmax] ) ) 
    
    cov_ens_obs = np.squeeze( cov_ens[:, obs_ii, obs_jj] )
    mod_lon_obs = np.squeeze( mod_lon[obs_ii, obs_jj] )
    mod_lat_obs = np.squeeze( mod_lat[obs_ii, obs_jj] )
    mod_z_obs   = np.squeeze( mod_z[obs_ii, obs_jj] ) 

    del cov_ens, mod_lon, mod_lat, mod_z, land_mask

    # Open process pool, send variables to GTM_OI_by_index, iterating over
    # elements of dom_index_list.

    pool = mp.Pool(n_proc)    
    print('Pool Open.')
    
    par_list = [(mod_lon_dom[ii], mod_lat_dom[ii], mod_z_dom[ii], 
                                 land_mask_dom[ii], mod_lon_obs, mod_lat_obs, mod_z_obs, 
                                 obs_z, Rvar, Brad, cov_ens_dom[ii], cov_ens_obs,
                                 ii) for ii in range(0, dL)]
    print('par_list created')
    par_out = pool.starmap(GTM_OI_by_index, par_list)
    print('Parallel routine finished.')
    pool.close(); print('Parallel pool closed.')
    pool.join(); print('Parallel pool joined.')
    
#   Output from parallel call will be in a list: [xa, inc, dom_index_list].
#   Take list and reconstruct full domain for xa and inc.
    xa_global = np.zeros((nr , nc)) * np.nan
    inc_global = np.zeros((nr , nc)) * np.nan
    #par_out = list(par_out)
    for ii in range(0,dL):
        out_tmp = par_out[ii]
        xa_tmp = out_tmp[0]; inc_tmp = out_tmp[1]; dom_index = out_tmp[2]
        ind_tup = dom_index_list[dom_index]
        xa_global[ind_tup[0]:ind_tup[1], ind_tup[2]:ind_tup[3]] = xa_tmp
        inc_global[ind_tup[0]:ind_tup[1], ind_tup[2]:ind_tup[3]] = inc_tmp
    
    return xa_global, inc_global

##############################
 ####GENERAL ROUTINES########
##############################
def gen_H_nn( xbL, obs_ind ):
    '''
    Generates a tangent linear matrix using nearest neighbour interpolation. 
    I.E. H(x) = y.
    #
 IN # xbL     :: Length of background state vector
 IN # obs_ind :: Indices of observation locations inside state vector
    '''
    
    obsL = len(obs_ind)
    H = np.zeros((obsL, xbL))
    H[np.arange(0,obsL), obs_ind] = 1
    
    return H

def gen_R( var ):
    '''
    Generates an observation covariance matrix assuming that observations are
    all independent (therefore only variances are required). Variances are
    placed on the diagonal of square matrix R.
    #
 IN # var :: vector of ordered observation error variances.
    '''
    
    R = np.diag(var)
    
    return R

def gen_B_para(D, c1, c2):
    '''
    Generates a parametric background covariance matrix
    '''
    
    B = c2*np.exp(-D**2/c1**2)
    
    return B

def OI(xb, y, B, H, R):
    '''
    Handles  the OI algorithm assuming all matrices are known and stored. 
    Generic routine.
    '''
    
    ##Check for NaN in observations and remove from y, H and R if necessary
    ynans = np.where(~np.isnan(y))
    y = np.array(y[ynans])
    H = np.squeeze(H[ynans,:])
    R = np.squeeze(R[ynans, :])
    R = np.squeeze(R[:,ynans])


    K = kalman_weights(B, H, R)
    innov = get_innov(xb, y, H)
    inc = np.matmul(K,innov)
    xa = xb + inc

    return xa, innov, inc

def kalman_weights(B, H, R):
    '''
    Calculates the kalman gain matrix, assuming all matrices are known.
    Generic routine.
    '''
    
    BH = np.matmul(B, np.transpose(H))
    HBH = np.matmul(H,BH)
    RHS = HBH + R
    RHSinv = np.linalg.inv(RHS)
    K = np.matmul(BH,RHSinv)
    
    return K

def OI2(xb, innov, BH, HBH, R):
    '''
    Handles  the OI algorithm assuming all matrices are known and stored. 
    Modified version of OI that assumes precalculated BH^T, HBH^T and 
    innovations. Useful for when observations lie outside of the assimilation
    domain.
    '''

    K = kalman_weights2(HBH,BH, R)
    inc = np.matmul(K,innov)
    xa = xb + inc

    return xa, inc

def kalman_weights2(HBH, BH, R):
    '''
    Calculates the kalman gain matrix, assuming all matrices are known.
    Modified version of OI that assumes precalculated BH^T, HBH^T and 
    innovations. Useful for when observations lie outside of the assimilation
    domain.
    '''
    RHS = HBH + R
    RHSinv = np.linalg.inv(RHS)
    K = np.matmul(BH,RHSinv)
    
    return K

def gen_A():
    
    return None

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
    
    D = dbg.dist_matrix_haversine(lon,lat,lon,lat)
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

def superobs_old(lon, lat, obs, crit_dist = 50, method = 'interp', obs_id=None):
    
    if obs_id is None:
        obs_id = np.zeros(lon.shape)
    
    D = dbg.dist_matrix_haversine(lon,lat,lon,lat)
    np.fill_diagonal(D,10000)
    crit = np.nanmin(D)
    
    lon_sup = np.array(lon)
    lat_sup = np.array(lat)
    obs_sup = np.array(obs)
    
    while crit < crit_dist:
        minii = np.unravel_index(np.argmin(D), D.shape)
        ii = minii[0]
        jj = minii[1]
        if method == 'interp':
            lon_sup[ii] = np.nanmean([lon_sup[ii], lon_sup[jj]])
            lat_sup[ii] = np.nanmean([lat_sup[ii], lat_sup[jj]])
            obs_sup[ii] = np.nanmean([obs_sup[ii], obs_sup[jj]])
            lon_sup = np.delete(lon_sup, jj)
            lat_sup = np.delete(lat_sup, jj)
            obs_sup = np.delete(obs_sup, jj)
        elif method == 'delete':
            lon_sup = np.delete(lon_sup, jj)
            lat_sup = np.delete(lat_sup, jj)
            obs_sup = np.delete(obs_sup, jj)
            obs_id = np.delete(obs_id, jj)
        
        D = dbg.dist_matrix_haversine(lon_sup, lat_sup, lon_sup, lat_sup)
        np.fill_diagonal(D,1000)
        crit = np.min(D)
    
    return lon_sup, lat_sup, obs_sup, obs_id

def cov_ssh(a, g):
    '''
    Creates a covariance matrix of model SSH for a specified harmonic
    constituent. Used to estimate background error covariance in the absence
    of observations.
    ''' 
    ssh = gen_ssh(a,g)
    C = np.cov(ssh)
    return C

def cov_innov(a_mod, g_mod, a_obs, g_obs, H):
    '''
    Estimates ssh error covariances between observation locations,
    '''
    a_mod = np.matmul(H, a_mod)
    g_mod = np.matmul(H, g_mod)
    
    ssh_mod = gen_ssh(a_mod, g_mod)
    ssh_obs = gen_ssh(a_obs, g_obs)
    
    err = ssh_mod - ssh_obs
    C = np.corrcoef(err)
    
    return C

def gen_ssh(a, g, L=2500):
    
    a = np.array(a)
    g = np.array(g)
    
    dsL = len(a)
    ssh = np.zeros((dsL,L))
    #Remove nans
    a[np.where(np.isnan(a))] = 0
    g[np.where(np.isnan(g))] = 0
    #Build ssh
    for ii in range(0,dsL):
        ssh[ii,:] = a[ii]*np.cos(np.arange(0,500,0.2) + np.radians(g[ii]))
    
    return ssh

def get_innov(xb, y, H):
    
    Hxb = np.matmul(H,xb)
    innov = y - Hxb
    
    return innov

def debias_linear(var, p):
    
    var_out = var - p[1]
    var_out = var_out/p[0]
    
    return var_out

############################################|
##########REGION/MASK DEFINITIONS###########|
############################################|

def domain_mask_by_index(nr, nc, domain_id):
    
    domain_mask = np.zeros((nr,nc))
    
    if domain_id == 0: #TEST BOX
        domain_mask[400:500, 950:1050] = 1
    
    return domain_mask    

def domain_mask_by_custom(nr, nc, rows, cols):
    
    domain_mask = np.zeros((nr,nc))
    domain_mask[rows[0]:rows[1], cols[0]:cols[1]] = 1
    
    return domain_mask    

def extract_dom_by_mask(mask):
    
    locs = np.where(mask)
    rmin = np.min(locs[0])
    rmax = np.max(locs[0])
    cmin = np.min(locs[1])
    cmax = np.max(locs[1])
    
    return rmin, rmax, cmin, cmax

def extract_obs_var_by_mask(obs_ii, obs_jj, mask):
    
    tmp = mask[obs_ii, obs_jj]
    ind = np.where(tmp==1)
    
    return ind   

def extract_obs_by_radius(ctr_lon, ctr_lat, radius, obs_lon, obs_lat):
    
    obsL = len(obs_lon)
    D = dbg.dist_haversine(np.ones(obsL)*ctr_lon, np.ones(obsL)*ctr_lat, 
                           obs_lon, obs_lat)
    obs_use = np.where(D<radius)
    
    return np.array(obs_use)

def replace_model_var_by_mask(var, mask):
    
    nr, nc = mask.shape
    var_out = np.zeros((nr,nc))*np.nan
    locs = np.where(mask)
    rmin = np.min(locs[0])
    rmax = np.max(locs[0])
    cmin = np.min(locs[1])
    cmax = np.max(locs[1])
    var_out[rmin:rmax,cmin:cmax] = var
    return var_out

############################################|
#############Handler Routines###############|
############################################|
    
def assimilate_principal_1var(xb, xb_lon, xb_lat, xb_mask, 
                           y, y_lon, y_lat, y_id, cov_ens, ass_rad = 6000):
    '''
    #
    # xb      - Background arrays (field, lon, lat)
    # y       - Observation vectors (vector, lon, lat)
    # Rvar    - Observation error variance
    # ass_rad - assimilation radius (km)
    '''

    print('Initialising..')
    #_____________________ORGANISE INPUT______________________________________#
    
    #Define observation input arrays
    y = np.array( y )
    y_lon = np.array( y_lon )
    y_lat = np.array( y_lat )
    
    # Define background arrays
    xb_lon = np.array( xb_lon )
    xb_lat = np.array( xb_lat )
    xb = np.array( xb )
    
    # Dataset sizes
    nr, nc = xb.shape
    
    #___________________ASSIMILATION__________________________________________#

    # Remove NAN from obs
    ynans = ~np.isnan(y)
    y  = np.array( y[ynans] )
    y_lon = np.array( y_lon[ynans] )
    y_lat = np.array( y_lat[ynans] )
    print('NaNs removed: ' + str( sum( ~ynans )) )
        
    print('Generating observation-model indices')
    # Obs indices after superobservation processing
    y_ii, y_jj = dbg.geo_nn_indices(xb_lon, xb_lat, y_lon, y_lat, 
                                    mask = xb_mask)
    
    # Define Rvar
    n_obs = len(y_lon)
    Rvar  = np.ones(n_obs) * np.nanmean(np.abs(y)) * 0.002 
    print(Rvar[0])
    Rvar  = Rvar**2
        
    # Calculate and remove any strong linear bias from model dataset
#    fit = np.polyfit(y, xb[y_ii, y_jj],1)
#    print(fit)
#    rs = dbg.r_squared_lin(y, xb[y_ii, y_jj], fit)
#    if rs > 0.8: xb_db = debias_linear(xb, fit)   # xb_debiased
#    else:        xb_db = xb
    xb_db = xb
    #Temporary covariance ensemble - demean
    cov_ens = np.array(cov_ens)
    cov_ens = cov_ens - np.nanmean(cov_ens, axis=0)
        
    # Do the assimilation
    print('Calling parallel assimilation routine')
    xa, inc = GTM_OI_global_par(xb_lon, xb_lat, xb_db, xb_mask, 
                                y_lon, y_lat, y, ass_rad, Rvar, [200,200], 
                                cov_ens = cov_ens, obs_ii = y_ii, obs_jj = y_jj,
                                n_proc = 15)
    
    #___________________POST-ASSIM__________________________________________#
    print('Assim done')
    inn = y - xb[y_ii, y_jj]
    ea  = y - xa[y_ii, y_jj]    
    y_info = (y_lon, y_lat, y, inn, y_ii, y_jj)

    return xa, inc, y_info, ea

    
def assimilate_multivar(xb, xb_lon, xb_lat, xb_mask, cov_ens, 
                        y, y_lon, y_lat,  ass_rad = 3000, boxsize = [150,150],
                        start_var = 'auto'):
    '''
    # For assimilation with an arbitrary number of variables.
    # Input variables are tuples (xb...y_lat1). Each element of a tuple is
    # data for a different variable. There must be as many background states
    # as variables. There can be fewer observed variables than actually. If
    # this is the case, the input observation tuples must still have the same
    # length as #var. Where observations are not available, set the tuple
    # element to [].
    '''

    print('Initialising..')
    #_____________________ORGANISE INPUT______________________________________#
    
    n_var = len(xb)
    
    # Remove NAN from obs
    for ii in range(0,n_var):
        if len(y[ii]) > 0:
            ynans = ~np.isnan(y[ii])
            y[ii] = np.array( y[ii][ynans] )
            y_lon[ii] = np.array( y_lon[ii][ynans] )
            y_lat[ii] = np.array( y_lat[ii][ynans] )
    
    # Dataset size
    nr, nc = xb[0].shape;
    
    if start_var == 'auto':
        rowbnds = 'auto'
    else:
        #Python indexing!!
        rowbnds = [start_var*nr, nr*n_var]
        print(rowbnds)
    
    #___________________PROCESSING___________________________________________#
        
    print('Generating observation-model indices')
    # Obs indices for each variable.
    y_ii = [[] for ii in range(0,n_var)]
    y_jj = [[] for ii in range(0,n_var)]
    for ii in range(0,n_var):
        if len(y_lon[ii]) > 0:
            tmpii, tmpjj = dbg.geo_nn_indices(xb_lon[ii], xb_lat[ii], 
                                              y_lon[ii], y_lat[ii], 
                                              mask = xb_mask[ii])
            y_ii[ii] = tmpii + ii*nr
            y_jj[ii] = tmpjj
        else:
            y_ii[ii] = []; y_jj[ii] = []
    
    print('Defining obs error variance')
    # Define Rvar
    n_obs = [len(ii) for ii in y_lon]
    Rvar = [[] for ii in range(0,n_var)]
    for ii in range(0,n_var):
        if len(y_lon[ii]) > 0:
            mn = np.nanmean(np.abs(y[ii]))
            Rvar[ii] = np.ones(n_obs[ii])*mn*0.002
            Rvar[ii] = Rvar[ii]**2
    
    # Calculate and remove any strong linear bias from model dataset
    #fit = np.polyfit(y, xb[y_ii, y_jj],1)
    #print('Fit: ' + str(fit))
    #rs = dbg.r_squared_lin(y, xb[y_ii, y_jj], fit)
    #if rs > 0.8: 
    #    print('Linear Adjustment, rs=' + str(rs))
    #    xb = debias_linear(xb, fit)   # xb_debiased
    #else:        print('Linear Adjustment, rs=' + str(rs))

    # Concatenate arrays for passing to assimilation routine
    xb = np.concatenate(xb, axis=0)
    xb_lon = np.concatenate(xb_lon, axis=0)
    xb_lat = np.concatenate(xb_lat, axis=0)
    xb_mask = np.concatenate(xb_mask, axis=0)
    y = np.concatenate(y)
    y_lon = np.concatenate(y_lon)
    y_lat = np.concatenate(y_lat)
    y_ii = np.concatenate(y_ii).astype(int)
    y_jj = np.concatenate(y_jj).astype(int)
    Rvar = np.concatenate(Rvar)
    cov_ens = np.concatenate(cov_ens, axis = 1)
    print(xb.shape)
        
    #___________________ASSIMILATION__________________________________________#
    # Do the assimilation
    print('Calling parallel assimilation routine')
    xa, inc = GTM_OI_global_par(xb_lon, xb_lat, xb, xb_mask, 
                                y_lon, y_lat, y, ass_rad, Rvar, boxsize, 
                                cov_ens = cov_ens, obs_ii = y_ii, obs_jj = y_jj,
                                n_proc = 20, rowbnds = rowbnds)
    
    #___________________POST-ASSIM__________________________________________#
#    print('Assim done, reshaping.')
    xa_out = [[] for ii in range(0,n_var)]
    inc_out = [[] for ii in range(0,n_var)]
    for ii in range(0,n_var):
        xa_out[ii] = xa[ii*nr:(ii+1)*nr]
        inc_out[ii] = inc[ii*nr:(ii+1)*nr]

    return xa_out, inc_out
    
    
def process_obs(y_lon, y_lat, y, obs_id, nemo_lon, nemo_lat, nemo_mask,
                fes_lon, fes_lat, fes_z, fes_mask, sobs_rad = 50,
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
        y_ii, y_jj = dbg.geo_nn_indices(nemo_lon, nemo_lat, y_lon, y_lat, 
                                        mask = nemo_mask)
    nlon = nemo_lon[y_ii, y_jj]; nlat = nemo_lat[y_ii, y_jj]
    dist = dbg.dist_haversine(nlon, nlat, y_lon, y_lat)
    nanii = dist>30
    y[nanii] = np.nan
    print('Points with bad model grid correspondence: ' + str(sum(nanii)))
    
    # Compare to FES data and take only obs that are within 1/2 error st. dev
    if fes_ii is None:
        fes_ii, fes_jj = dbg.geo_nn_indices(fes_lon, fes_lat, y_lon, y_lat, 
                                    mask = fes_mask)
    y_fes = fes_z[fes_ii, fes_jj]
    ef = y - y_fes
    estd = np.nanstd(ef); emean = np.nanmean(ef)
    aef = np.abs(ef - emean)
    nanii = aef > 0.5*estd
    y[nanii] = np.nan
    print('Comparison to FES data made: ' + str( sum(nanii )) + ' removed.')
    
    # Superobs
    y_lon, y_lat, y, obs_id = superobs(y_lon, y_lat, y, crit_dist= sobs_rad, 
                               method='delete', obs_id = obs_id)
    print('Superobs processed.')
    
    return y_lon, y_lat, y, obs_id
    
    
