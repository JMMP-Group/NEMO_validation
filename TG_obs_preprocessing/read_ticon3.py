"""
File to preprocess the observational harmonic data
TICON-3

Hart-Davis, Michael G; Dettmering, Denise; Seitz, Florian (2022): TICON-3: Tidal Constants based on GESLA-3 sea-level records from globally distributed tide gauges including gauge type information (data) [dataset]. PANGAEA, https://doi.org/10.1594/PANGAEA.951610



    # # # # # READ MDP HARMONICS # # # # #
    tmp = dbr.read_mdp_tg_harm(fn_mdp_data, mdp_data_types, cc)
    mdp_a = tmp[0]; mdp_g = tmp[1]; mdp_lon = tmp[2]+360; mdp_lat = tmp[3]
    mdp_z1 = tmp[4]; mdp_z2 = tmp[5]
    mdp_lon[mdp_lon>433] = mdp_lon[mdp_lon>433] - 360

    const = ['SA','SSA','MF','MM','MSF','Q1','O1','P1','S1','K1','J1','2N2','MU2',
         'N2','NU2','M2','MKS2','L2','T2','S2','R2','K2','M3','MN4','M4','MS4',
         'S4','M6','M8']
    n_const = len(const)
    cc = const[ii]

    
    Consituents in TICON-3:
    array(['M2', 'K1', 'N2', 'O1', 'P1', 'Q1', 'K2', 'S2', 'S1', 'SA', 'T2',
       'MF', 'MM', '2N2', 'M4', 'J1', 'SSA', 'MSF', 'MSQ', 'EP2', 'L2',
       'M3', 'R2 ', 'MI2', 'MTM', 'NI2', 'LM2', 'MN4', 'MS4', 'MKS', 'N4',
       'M6', 'M8', 'S4', '2Q1', 'OO1', 'S3', 'MA2', 'MB2', 'M1'],
      dtype=object)

    TICON-3 is "missing": ['R2', 'MU2', 'MKS2', 'NU2']
    """

# useful packages
import numpy as np
import pandas as pd

fn_ticon3_data = "/Users/jelt/Downloads/TICON_3/TICON_3.txt" # path to the TICON-3 data
cc = 'M2' # constituent to select, e.g. 'M2'

def read_ticon3(fn_ticon3_data, cc, gauge_type='Coastal'):
    """
    Reads the TICON-3 data file and returns the amplitude, phase, longitude, latitude,
     z1 and z2 components for a specified constituent and gauge type.

    If amplitude/phase are flagged or not present, they will be set to NaN.
    Latitude and longitude will be set to zero for this request.

    jelt 2025-06-4
    """
    
    
    # load data
    ptd = "/Users/jelt/Downloads/TICON_3/" # path to the TICON-3 data
    fnm ="TICON_3.txt" # TICON-3 file
    df = pd.read_csv(fn_ticon3_data, sep='\t', header=None)
    
    # initialize variables
    lat = 0; lon = 0; amp = np.nan; pha = np.nan; z1 = np.nan; z2 = np.nan

    try:
        # select a constituent
        indx_cons = np.where(df[2] == cc)
        df = df.iloc[indx_cons]
        # select tide gauge type [options are: River, Lake and Coastal]
        indx_type = np.where(df[13] == gauge_type)
        df = df.iloc[indx_type]
        lon = np.array(df[1]) # assign longitude
        lat = np.array(df[0]) # assign latitude
        amp = np.array(df[3]) # assign amplitude
        pha = np.array(df[4]) # assign phase
        con = np.array(df[2]) # assign consituent name
        # Also change interval phase lies in from 0 -> 360 to -180 -> 180.
        pha[pha>180] = pha[pha>180] - 360
        # Convert amplitudes and phases to z1 and z2
        z1, z2 = dbg.polar2cart(amp, pha, degrees=True)
    except Exception as e:
        pass
    return amp, pha, lon, lat, z1, z2

tmp = read_ticon3(fn_ticon3_data, cc)

tidcon_amp = tmp[0]  # amplitude
tidcon_pha = tmp[1]  # phase
tidcon_lon = tmp[2]  # longitude
tidcon_lat = tmp[3]  # latitude
tidcon_z1 = tmp[4]  # z1 component
tidcon_z2 = tmp[5]  # z2 component