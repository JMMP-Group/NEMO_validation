import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/work/dbyrne/code/COAsT')
import coast
from validate_transects import extract_transects_t


fn_nemo_domain = "/projectsa/NWSrivers/AMM7/INPUTS/mesh_mask.nc"
fn_nemo_data = "/projectsa/NWSrivers/AMM7/REF/*daily.grid_T*"
dn_out = "/projectsa/NWSrivers/AMM7/analysis/REF/transects"


lat0_B = np.array( [58, 60, 60, 59, 59, 59, 58, 57, 54, 54, 52, 51, 50, 50, 52, 53, 53, 55] )
lat0_s = np.array( [32, 41, 41, 55, 17, 17, 10, 8, 15, 15, 25, 19, 49, 25, 17, 27, 15, 0])

lon0_B = np.array( [3, 0, 3, 1, 2, 3, 6, 2, 0, 4, 1, 1, 1, 4, 6, 6, 6, 6])
lon0_s = np.array( [10, 55, 0, 15, 30, 20, 32, 15, 35, 14, 40, 20, 25, 25, 10 ,20, 20, 20])

lat1_B = np.array( [59, 60, 60, 59, 59, 59, 57, 57, 54, 53, 52, 50, 49, 48, 52, 53, 55, 55])
lat1_s = np.array( [17, 41, 41, 17, 17, 17, 8, 8, 15, 15, 25 ,55, 45, 35, 17, 27, 0, 0])

lon1_B = np.array( [2, 3, 6, 2, 3, 6, 5, 5, 4, 5, 4, 1, 1, 4, 5, 4, 6, 8])
lon1_s = np.array( [30, 0, 0, 30, 20, 30, 0, 0, 14, 30, 35, 42, 25, 25, 0, 35, 20, 40])


A_lat = lat0_B + lat0_s/60
A_lon = lon0_B + lon0_s/60
B_lat = lat1_B + lat1_s/60
B_lon = lon0_B + lon0_s/60

transect_names = ['Pentland_Firth', 'Shetland_N', 'Sognesjoen', 'Shetland_S',
                  'Orkney', 'Utsira', 'Lista', 'Aberdeen', 'Flamborough_Hd',
                  'Terschelling','Noordwijk','Dover_Strait','Cherbourgh',
                  'Plymouth','Rosslare','Dublin_Holyhead','Rottumerplaat', 'Sylt']

extract_transects_t(fn_nemo_data, fn_nemo_domain, A_lat, A_lon, 
                    B_lat, B_lon, dn_out, transect_names=transect_names)