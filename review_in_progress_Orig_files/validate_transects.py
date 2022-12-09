# UNCOMMENT IF USING A DEVELOPMENT VERSION OF COAST
import sys
sys.path.append('/Users/dbyrne/code/COAsT/')

import coast
import coast.general_utils as coastgu
import coast.plot_util as pu
import numpy as np
import xarray as xr
import os.path
from dask.diagnostics import ProgressBar

def write_ds_to_file(ds, fn, **kwargs):
    ''' 
    Simple netcdf writing routine which checks if file already exists first 
    '''
    if os.path.exists(fn):
        os.remove(fn)
    with ProgressBar():
        ds.to_netcdf(fn, **kwargs)

def extract_transects_t(fn_nemo_data_t, fn_nemo_domain,
                      A_lon, A_lat, B_lon, B_lat,
                      dn_out, fn_nemo_data_u=None, fn_nemo_data_v=None,
                      fn_out_pre = 'transect', transect_names=[],
                      nemo_chunks = {'time_counter':250}, multiple=True):
    
    nemo = coast.NEMO(fn_nemo_data_t, fn_nemo_domain, 
                      grid_ref='t-grid', chunks=nemo_chunks, multiple=multiple)
    
    if fn_nemo_data_u is not None:
        nemo_u = coast.NEMO(fn_nemo_data_u, fn_nemo_domain, 
                      grid_ref='u-grid', chunks=nemo_chunks, multiple=multiple)
        nemo_v = coast.NEMO(fn_nemo_data_v, fn_nemo_domain, 
                      grid_ref='v-grid', chunks=nemo_chunks, multiple=multiple)
        nemo_f = coast.NEMO(fn_domain = fn_nemo_domain, grid_ref='f-grid', chunks={})
    
    n_transect = len(A_lon)
    
    for ii in range(n_transect):
        print('Transect {0}/{1}'.format(ii, n_transect))
        tran = coast.Transect_t( nemo, (A_lon[ii], A_lat[ii]), (B_lon[ii], B_lat[ii]))
        
        if fn_nemo_data_u is not None:
            tran_f = coast.Transect_f(nemo_f, (A_lon[ii], A_lat[ii]), (B_lon[ii], B_lat[ii])) 
            tran_f.calc_flow_across_transect(nemo_u,nemo_v)
            tran = xr.merge((tran, tran_f))
        
        ds = tran.data
        ds['x_ind'] = (['r_dim'], tran.x_ind)
        ds['y_ind'] = (['r_dim'], tran.y_ind)
        ds['A_lon'] = A_lon[ii]
        ds['A_lat'] = A_lon[ii]
        ds['B_lon'] = B_lon[ii]
        ds['B_lat'] = B_lon[ii]
        
        if len(transect_names) == 0:
            fn_out = '{0}_{1}N{2}E_{3}N{4}E.nc'.format(fn_out_pre, A_lat[ii], A_lon[ii], B_lat[ii], B_lon[ii])
        else:
            fn_out = '{0}_{1}.nc'.format(fn_out_pre, transect_names[ii])
        fn_out = os.path.join(dn_out, fn_out)
        write_ds_to_file(ds, fn_out)


    



