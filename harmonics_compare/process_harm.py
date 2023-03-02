import os
os.environ['DN_OUT'] ="/Users/jelt/Downloads/SENEMO/"
os.environ['COAST_REPO'] = "/home/users/jelt/GitHub/COAsT"

from config import config
import sys

config = config() # initialise variables in python

# IF USING A DEVELOPMENT BRANCH OF COAST, ADD THE REPOSITORY TO PATH:
sys.path.append(config.coast_repo)

import coast
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
#import datetime
import pandas as pd
import os
from coast import crps_util as cu
import numpy as np

import time
#from validate_ssh_tg_hourly import extract_ssh, analyse_ssh, plot_single_cfg
import numpy as np


def modulo_align_phases_to_x(x, y, deg=True):
    """ Use modulo arithmetic to align y values with x in order to be able to calculate regression """
    if deg == True:
        circ = 360
        x = (x + 180) % 360 - 180
    else:
        circ = 2 * np.pi
        x = (x * 180 / np.pi + 180) % 360 - 180
        x = x * np.pi / 180.

    I_y_large = y - x > 0.5 * circ
    I_y_small = x - y > 0.5 * circ

    y[I_y_large] = y[I_y_large] - circ
    y[I_y_small] = y[I_y_small] + circ

    return x, y

def pearson_correl_coef(z1obs, z2obs, z1mod, z2mod):
    """
    sum( <zobs_i, zmod_i> )/ ( rms(zobs) rms(zmod) )
        for vectors zobs_i = (z1obs_i, z2obs_i))
        and zero mean zobs and zmod
    """
    inner_product = np.nanmean( z1obs*z1mod + z2obs*z2mod )
    return inner_product / (rms_abs_amp(z1obs,z2obs) * rms_abs_amp(z1mod,z2mod))

def rms_abs_amp(z1,z2):
    """
    root mean square vector amplitude
    z1=a.cos(theta), z2=a.sin(theta)
    """
    return np.sqrt(np.nanmean( z1**2 + z2**2 ))

def rms_abs_error(z1obs, z2obs, z1mod, z2mod):
    """
    root mean square absolute error = sqrt( sum |zmod - zobs|^2 )
    """
    return rms_abs_amp( z1obs - z1mod, z2obs - z2mod )

def r_squared_lin(x, y, fit):
    '''For calculating r-squared of a linear fit. Fit should be a python polyfit
    object'''

    fity = fit(x)
    diff = (y - fity) ** 2
    ybar = np.nanmean(y)
    ymybar = (y - ybar) ** 2

    SStot = np.nansum(ymybar)
    SSres = np.nansum(diff)

    R2 = 1 - SSres / SStot

    return R2


def rmsd(x1, x2):
    # Calculated root mean square difference of 2 given vectors of values x1 and x2. One is likely obs
    # and the other model output.

    diff = x1 - x2
    diffs = diff ** 2
    msd = np.nanmean(diffs)
    rmsd = np.sqrt(msd)

    return rmsd


def score(y, x):
    """ return slope, r2 score, and rms difference """
    I = ~np.isnan(x) * ~np.isnan(y)

    reg = np.polyfit(x[I], y[I], 1)
    print(f"regression model: {reg}")
    slope = reg[0]

    predict = np.poly1d(reg)
    score = r_squared_lin(x[I], y[I], predict)
    rms = rmsd(x[I], y[I])
    # print(f"r2 score: {score}")
    return float(format(slope, '.2f')), float(format(score, '.2f')), float(format(rms, '.2f'))


def plot_scatter_score(X1, Y1, X2, Y2, subtitle_str=["", ""], title_str: str = "", yex=True):
    y_d = Y1[ind_deep]
    x_d = X1[ind_deep]
    sl_d, r2_d, rms_d = score(y_d, x_d)

    y_s = Y1[~ind_deep]
    x_s = X1[~ind_deep]
    sl_s, r2_s, rms_s = score(y_s, x_s)

    plt.close('all')
    plt.figure()
    plt.subplot(2, 1, 1)
    plt.plot(x_d, y_d, 'k.', label=f">200m (m:{sl_d}, r2:{r2_d}, rms:{rms_d})")
    plt.plot(x_s, y_s, 'g+', label=f"<200m (m:{sl_s}, r2:{r2_s}, rms:{rms_s})")
    if yex:
        [xmin, xmax] = plt.gca().get_xlim()
        [ymin, ymax] = plt.gca().get_ylim()
        plt.plot([max(xmin, ymin), min(xmax, ymax)], [max(xmin, ymin), min(xmax, ymax)], 'k--')
        # plt.plot([-1,1],[-1,1],'k--')
    plt.xlabel('obs')
    plt.ylabel('mod')
    plt.title(subtitle_str[0])
    plt.legend()

    y_d = Y2[ind_deep]
    x_d = X2[ind_deep]
    sl_d, r2_d, rms_d = score(y_d, x_d)

    y_s = Y2[~ind_deep]
    x_s = X2[~ind_deep]
    sl_s, r2_s, rms_s = score(y_s, x_s)

    plt.subplot(2, 1, 2)
    plt.plot(x_d, y_d, 'k.', label=f">200m (m:{sl_d}, r2:{r2_d}, rms:{rms_d})")
    plt.plot(x_s, y_s, 'g+', label=f"<200m (m:{sl_s}, r2:{r2_s}, rms:{rms_s})")
    if yex:
        [xmin, xmax] = plt.gca().get_xlim()
        [ymin, ymax] = plt.gca().get_ylim()
        plt.plot([max(xmin, ymin), min(xmax, ymax)], [max(xmin, ymin), min(xmax, ymax)], 'k--')
        # plt.plot([-1,1],[-1,1],'k--')
    plt.xlabel('obs')
    plt.ylabel('mod')
    plt.title(subtitle_str[1])
    plt.legend()
    plt.suptitle(title_str)
    plt.tight_layout()
    #plt.show()
    plt.savefig(config.dn_out+"PROCESSED/FIGS/M2_fit_" + title_str + ".png")


# Load data as tidegauge objects
## Harmonise definitions: negate M2y (and phase) in NEMO
################################
obs = coast.Tidegauge(dataset=xr.open_dataset(config.dn_out+"PROCESSED/obs_extracted.nc"))
obs.dataset['longitude'] = ((obs.dataset.longitude%360) + 180)%360 - 180  # set longitude : [-180,180]
obs.dataset['G'] = (obs.dataset.G+180)%360-180  # set phase: -180,180

fes = coast.Tidegauge(dataset=xr.open_dataset(config.dn_out+"PROCESSED/FES2014_extracted.nc"))
fes.dataset['G'] = (fes.dataset.G+180)%360-180  # set phase: -180,180

gs1p1 = coast.Tidegauge(dataset=xr.open_dataset(config.dn_out+"PROCESSED/GS1p1_tide_extracted.nc"))
gs1p1.dataset['G'] = -(gs1p1.dataset.G+180)%360-180  # set phase: -180,180
gs1p1.dataset['M2y'] = - gs1p1.dataset.M2y
gs1p1.dataset = gs1p1.dataset.drop_dims(["nvertex", "z_dim"])  # drop unwanted dimensions and associated variables

gs1p2 = coast.Tidegauge(dataset=xr.open_dataset(config.dn_out+"PROCESSED/GS1p2_full_extracted.nc"))
gs1p2.dataset['G'] = -(gs1p2.dataset.G+180)%360-180  # set phase: -180,180
gs1p2.dataset['M2y'] = - gs1p2.dataset.M2y
gs1p2.dataset = gs1p2.dataset.drop_dims(["nvertex", "z_dim"])  # drop unwanted dimensions and associated variables



# Compare maps
obs.plot_on_map_multiple([obs], color_var_str="G")
#plt.savefig(config.dn_out+"PROCESSED/FIGS/obs_phase_on_map.png")

if(0):
    # definitions for M2x and M2y are opposite to NEMO, but work:
    fig, [ax0, ax1] = plt.subplots(ncols=2)
    fes.dataset.M2x.plot(ax=ax0)
    gs1p1.dataset.M2x.plot(ax=ax1)
    plt.show()


# Phase alignment
ylims = [-80,80] #[45,60]
xlims = [-180,180] #[-60, -25]

plt.close('all')
plt.figure()
plt.subplot(2,2,1)
plt.scatter(obs.dataset.longitude, obs.dataset.latitude, c=obs.dataset['G'], s=10 )
plt.colorbar()
plt.title('obs')
plt.xlim(xlims)
plt.ylim(ylims)

plt.subplot(2,2,2)
plt.scatter(fes.dataset.longitude, fes.dataset.latitude, c=fes.dataset['G'], s=10 )
plt.colorbar()
plt.title('FES')
plt.xlim(xlims)
plt.ylim(ylims)

plt.subplot(2,2,3)
plt.scatter(gs1p1.dataset.longitude, gs1p1.dataset.latitude, c=gs1p1.dataset['G'], s=10 )
plt.colorbar()
plt.title('GS1P1')
plt.xlim(xlims)
plt.ylim(ylims)

plt.subplot(2,2,4)
plt.scatter(gs1p2.dataset.longitude, gs1p2.dataset.latitude, c=gs1p2.dataset['G'], s=10 )
plt.colorbar()
plt.title('GS1P2')
plt.xlim(xlims)
plt.ylim(ylims)
#plt.show()
plt.savefig(config.dn_out+"PROCESSED/FIGS/scatter_phase_on_map.png")


#%% Align datasets
# TidegaugeAnalysis routine would work well if only _two_ dataarray had to be aligned:
# tganalysis = coast.TidegaugeAnalysis()
# for var in ['M2x', 'M2y']:
#    obs_new, mod_new = tganalysis.match_missing_values(tg_obs.dataset[var], tg_fes.dataset[var])
#    tg_obs.dataset[var] = obs_new.dataset[var]
#    tg_fes.dataset[var] = mod_new.dataset[var]
#
# But instead we do it by hand


ind1 = np.isnan(obs.dataset.M2x.values)
ind2 = np.isnan(obs.dataset.M2y.values)
ind3 = np.isnan(fes.dataset.M2x.values)
ind4 = np.isnan(fes.dataset.M2y.values)
ind5 = np.isnan(gs1p1.dataset.M2x.values)
ind6 = np.isnan(gs1p1.dataset.M2y.values)
ind7 = np.isnan(gs1p2.dataset.M2x.values)
ind8 = np.isnan(gs1p2.dataset.M2y.values)

I = ~ind1 * ~ind2 * ~ind3 * ~ind4 * ~ind5 * ~ind6 * ~ind7 * ~ind8

#obs.dataset.M2x[~I] = np.nan
#obs.dataset.M2y[~I] = np.nan
#fes.dataset.M2x[~I] = np.nan
#fes.dataset.M2y[~I] = np.nan
#gs1p1.dataset.M2x[~I] = np.nan
#gs1p1.dataset.M2y[~I] = np.nan
#gs1p2.dataset.M2x[~I] = np.nan
#gs1p2.dataset.M2y[~I] = np.nan
obs.dataset = obs.dataset.where(I)
fes.dataset = fes.dataset.where(I)
gs1p1.dataset = gs1p1.dataset.where(I)
gs1p2.dataset = gs1p2.dataset.where(I)

print(ind1.flatten().sum())
print(ind2.flatten().sum())
print(ind3.flatten().sum())
print(ind4.flatten().sum())
print(ind5.flatten().sum())
print(ind6.flatten().sum())
print(ind7.flatten().sum())
print(ind8.flatten().sum())
print((~I).flatten().sum())


#%%  Plot distributions of depth at observation locations
plt.close('all')
plt.figure()
plt.plot( np.sort(np.log10(gs1p1.dataset.bathymetry.values)) )
#plt.plot([0,500],[np.log10(2000), np.log10(2000)], 'm', label="2000m")
plt.plot([0,500],[np.log10(200), np.log10(200)], 'r', label="200m")
#plt.plot([0,500],[np.log10(50), np.log10(50)], 'g', label="50m")
#plt.plot([0,500],[np.log10(15), np.log10(15)], 'k', label="15m")
plt.xlabel('count')
plt.ylabel('log10(depth)')
plt.title("distribution of depths at observation sites")
plt.legend()
plt.savefig(config.dn_out+"PROCESSED/FIGS/dist_bathy.png")



#%% ## Plot amplitude and phase errors (deg)

for EXP in ["GS1P1", "GS1P2", "FES2014"]:

    if EXP == "GS1P1": tg_mod = gs1p1
    if EXP == "GS1P2": tg_mod = gs1p2
    if EXP == "FES2014": tg_mod = fes

    # separate observations by depth
    ind_deep = tg_mod.dataset.bathymetry.values > 200

    X1, Y1 = obs.dataset.M2x, tg_mod.dataset.M2x
    X2, Y2 = obs.dataset.M2y, tg_mod.dataset.M2y

    plot_scatter_score(X1,Y1, X2,Y2, subtitle_str=["M2x (m)","M2y (m)"], title_str=EXP+"_complex")


#%% ## Plot amplitude and phase errors (deg)

for EXP in ["GS1P1", "GS1P2", "FES2014"]:

    if EXP == "GS1P1": tg_mod = gs1p1
    if EXP == "GS1P2": tg_mod = gs1p2
    if EXP == "FES2014": tg_mod = fes

    # separate observations by depth
    ind_deep = tg_mod.dataset.bathymetry.values > 200

    X1, Y1 = obs.dataset.A, tg_mod.dataset.A
    X2, Y2 = modulo_align_phases_to_x(obs.dataset.G, tg_mod.dataset.G, deg=True)  # preprocess w/ modulo arithmetic

    plot_scatter_score(X1,Y1, X2,Y2, subtitle_str=["Amplitude (m)","Phase (deg)"], title_str=EXP+"_amp_pha")


#%% Compute Taylor diagram stats
count = 0
for subset in ['shal', 'deep']:
    try: del II, z1obs, z2obs, z1mod, z2mod
    except: pass
    if subset == 'deep':
        # separate observations by depth
        II = fes.dataset.bathymetry.values > 200
    elif subset == 'shal':
        II = fes.dataset.bathymetry.values <= 200
    else:
        print(f"Not expecting that {subset}")

    z1obs, z2obs = obs.dataset.M2x[II], obs.dataset.M2y[II]
    if(0):
        # %%  Plot distributions of depth at observation locations
        plt.close('all')
        plt.figure()
        plt.plot(np.sort(obs.dataset.A[I].values))
        plt.xlabel('count')
        plt.ylabel('M2 Amp (m))')
        plt.title("distribution of A at observation sites")
        plt.legend()
        plt.savefig(config.dn_out + "PROCESSED/FIGS/dist_obsA_"+subset+".png")


    # Obs
    if count == 0:
        R = np.array([1])  # R for obs
        rms_amp = np.array([rms_abs_amp(z1obs, z2obs)])
        rms_err = np.array([0])  # err for obs
    else:
        R = np.hstack((R, 1))
        rms_amp = np.hstack((rms_amp, rms_abs_amp(z1obs, z2obs)))
        rms_err = np.hstack((rms_err, 0))

    # FES
    z1mod, z2mod = fes.dataset.M2x[II], fes.dataset.M2y[II]
    # obs_new, mod_new = tganalysis.match_missing_values(obs.dataset.M2x, tg_fes.dataset.M2x)
    # z1obs_new = obs_new.dataset.M2x
    # z1mod_new = mod_new.dataset.M2x
    # obs_new, mod_new = tganalysis.match_missing_values(obs.dataset.M2y, tg_fes.dataset.M2y)
    # z2obs_new = obs_new.dataset.M2y
    # z2mod_new = mod_new.dataset.M2y

    R = np.hstack((R, pearson_correl_coef(z1obs, z2obs, z1mod, z2mod)))
    rms_amp = np.hstack((rms_amp, rms_abs_amp(z1mod, z2mod)))
    rms_err = np.hstack((rms_err, rms_abs_error(z1obs, z2obs, z1mod, z2mod)))

    # GS1P1
    del z1mod, z2mod
    z1mod, z2mod = gs1p1.dataset.M2x[II], gs1p1.dataset.M2y[II]
    # z1obs_new, z1mod_new = tganalysis.match_missing_values(obs.dataset.M2x, tg_mx2.dataset.M2x)
    # z2obs_new, z2mod_new = tganalysis.match_missing_values(obs.dataset.M2y, tg_mx2.dataset.M2y)

    R = np.hstack((R, pearson_correl_coef(z1obs, z2obs, z1mod, z2mod)))
    rms_amp = np.hstack((rms_amp, rms_abs_amp(z1mod, z2mod)))
    rms_err = np.hstack((rms_err, rms_abs_error(z1obs, z2obs, z1mod, z2mod)))

    # GS1P2
    del z1mod, z2mod
    z1mod, z2mod = gs1p2.dataset.M2x[II], gs1p2.dataset.M2y[II]
    # z1obs_new, z1mod_new = tganalysis.match_missing_values(obs.dataset.M2x, tg_mv2.dataset.M2x)
    # z2obs_new, z2mod_new = tganalysis.match_missing_values(obs.dataset.M2y, tg_mv2.dataset.M2y)
    R = np.hstack((R, pearson_correl_coef(z1obs, z2obs, z1mod, z2mod)))
    rms_amp = np.hstack((rms_amp, rms_abs_amp(z1mod, z2mod)))
    rms_err = np.hstack((rms_err, rms_abs_error(z1obs, z2obs, z1mod, z2mod)))


    count = count + 1
label = ['obs:s', 'fes:s', 'gs1p1:s', 'gs1p2:s',
         'obs:d', 'fes:d', 'gs1p1:d', 'gs1p2:d']


print(f"R= {[format(R[i],'.2f') for i in range(len(R))]}")
print(f"rms_amp= {[format(rms_amp[i],'.2f') for i in range(len(R))]}")
print(f"rms_err= {[format(rms_err[i],  '.2f') for i in range(len(R))]}")
print(label)


## Check cosine rule consistency
A = rms_amp
B = rms_amp[0]
C = rms_err
costheta = R

for i in range(len(R)):
    print(f"{label[i]}: sqrt(A^2+B^2-2ABcos(theta))={np.sqrt(A[i]**2 + B**2 - 2*A[i]*B*costheta[i])}. C={C[i]}")
