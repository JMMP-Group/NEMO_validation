"""
This script has been developed using local data on my macbook.
The idea would be to port it to JASMIN, or somewhere once all the debugging is done.

Useage:
python process_harm.py
"""


import os
os.environ['DN_OUT'] ="/Users/jelt/Downloads/SENEMO/"
#os.environ['COAST_REPO'] = "/home/users/jelt/GitHub/COAsT"
os.environ['COAST_REPO'] = "/Users/jelt/GitHub/COAsT"


from config import config
import sys

config = config() # initialise variables in python

# IF USING A DEVELOPMENT BRANCH OF COAST, ADD THE REPOSITORY TO PATH:
sys.path.append(config.coast_repo)

import coast
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib
import cartopy.crs as ccrs
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
from taylor_harmonic_tide_plot import TaylorTide

constit_list = ['M2','S2','N2','K1','O1','Q1','P1','K2'] #,'M4','MS4']

def align_datasets(tg_obs, tg_mod):
    """

    Match the nans between two datasets
    """
    for constit in constit_list:
        ind1 = np.isnan(tg_obs.dataset[constit+"x"].values)
        ind2 = np.isnan(tg_obs.dataset[constit+"y"].values)
        ind3 = np.isnan(tg_mod.dataset[constit+"x"].values)
        ind4 = np.isnan(tg_mod.dataset[constit+"y"].values)
        I = ~ind1 * ~ind2 * ~ind3 * ~ind4

        tg_obs.dataset[constit+"x"][~I] = np.nan
        tg_obs.dataset[constit+"y"][~I] = np.nan
        tg_mod.dataset[constit+"x"][~I] = np.nan
        tg_mod.dataset[constit+"y"][~I] = np.nan

        #tg_obs.dataset = tg_obs.dataset.where(I)
        #tg_mod.dataset = tg_mod.dataset.where(I)
    return tg_obs, tg_mod



def amp_pha_from_re_im(creal, cimag):
    """
    Example usage: amp,pha = amp_pha_from_re_im(ds.M2x, ds.M2y)
    """
    cc=creal+cimag*1j
    if type(creal) == xr.DataArray:
        amp = xr.zeros_like(creal)
        pha = xr.zeros_like(creal)
        amp=amp.where( np.isnan(amp),  np.abs(cc))
        pha=pha.where( np.isnan(pha),  np.angle(cc)*180/np.pi)
    else:
        amp = np.abs(cc)
        pha = np.angle(cc)*180/np.pi
    return amp, pha

def re_im_from_amp_pha(amp, pha):
    """
    Assumes phase in degrees
    Example usage: amp,pha = amp_pha_from_re_im(ds.amp, ds.pha)
    """
    if np.max(np.abs(pha)) <= 2*np.pi:
        print(f"Warning. Check phase units. Expected degrees. Max: {np.max(pha)}. Min: {np.min(pha)}")
    re = amp*np.cos(pha*np.pi/180)
    im = amp*np.sin(pha*np.pi/180)
    return re, im

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


def plot_scatter_score(X1, Y1, X2, Y2, ind_deep, constit_str=None, subtitle_str=["", ""], title_str: str = "", yex=True):

    if constit_str==None:
        print(f"plot_scatter_score(): No constituent specified")
        return

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
    plt.savefig(config.dn_out+"PROCESSED/FIGS/"+constit_str+"_fit_" + title_str + ".png")


def plot_scatter_line(X1, Y1, X2, Y2, ind_deep, constit_str=None, subtitle_str=["", ""], title_str: str = "", yex=True):
    # scatter plot the vector differences between two sets of vectors
    if constit_str==None:
        print(f"plot_scatter_line(): No constituent specified")
        return

    mod_x_d = Y1[ind_deep]
    obs_x_d = X1[ind_deep]

    mod_x_s = Y1[~ind_deep]
    obs_x_s = X1[~ind_deep]

    mod_y_d = Y2[ind_deep]
    obs_y_d = X2[ind_deep]

    mod_y_s = Y2[~ind_deep]
    obs_y_s = X2[~ind_deep]


    plt.close('all')
    plt.figure()
    plt.subplot(1, 1, 1)
    plt.plot([obs_x_d, mod_x_d], [obs_y_d, mod_y_d], 'k-')#, label=f">200m")
    plt.plot([obs_x_s, mod_x_s], [obs_y_s, mod_y_s], 'g-')#, label=f"<200m")
    if yex:
        [xmin, xmax] = plt.gca().get_xlim()
        [ymin, ymax] = plt.gca().get_ylim()
        plt.plot([0, 0], [ymin, ymax], 'k--')
        plt.plot([xmin, xmax], [0, 0], 'k--')
    plt.xlabel('A.cos(theta)')
    plt.ylabel('A.sin(theta)')
    plt.title(subtitle_str[0])

    colors = ['black', 'green']
    lines = [Line2D([0], [0], color=c, linewidth=3, linestyle='solid') for c in colors]
    labels = [">200m", "<200m"]
    plt.legend(lines, labels)

    plt.suptitle(title_str)
    plt.tight_layout()
    #plt.show()
    plt.savefig(config.dn_out+"PROCESSED/FIGS/"+constit_str+"_fit_" + title_str + ".png")


def plot_all_complex_harmonic_stick_errors():
    # Loop over model runs and harmonic species.
    # Plot vector errors between complex model and obs harmonics: (z1,z2)_mod - (z1,z2)_obs
    for EXP in ["ZPS_TIDE", "MES_TIDE", "FES2014"]:
        for constit in constit_list:
            if EXP == "ZPS_TIDE": tg_mod = nemo1
            if EXP == "MES_TIDE": tg_mod = nemo2
            if EXP == "FES2014": tg_mod = fes

            # separate observations by depth
            try:
                ind_deep = tg_mod.dataset.bathymetry.values > 200
            except:
                ind_deep = nemo1.dataset.bathymetry.values > 200
                print(f"Issue with bathy in {EXP}. Use bathy from GS1p1")

            try:
                X1, Y1 = obs.dataset[constit + "x"], tg_mod.dataset[constit + "x"]
                X2, Y2 = obs.dataset[constit + "y"], tg_mod.dataset[constit + "y"]

                plot_scatter_line(X1, Y1, X2, Y2, ind_deep, constit_str=constit,
                                   subtitle_str=['delta: ' + constit], title_str=EXP + "_complex_stick")
            except:
                plt.close('all')
                plt.figure()
                plt.plot(0, 0);
                plt.title(EXP + ":" + constit)
                plt.savefig(config.dn_out + "PROCESSED/FIGS/" + constit + "_fail_" + EXP + "_complex_stick.png")
                print(f"Issue with {constit} in {EXP}")

def plot_all_complex_harmonic_errors():
    # Loop over model runs and harmonic species.
    # Plot z1 & z2 with best fit stats against obs
    for EXP in ["ZPS_TIDE", "MES_TIDE", "FES2014"]:
        for constit in constit_list:
            if EXP == "ZPS_TIDE": tg_mod = nemo1
            if EXP == "MES_TIDE": tg_mod = nemo2
            if EXP == "FES2014": tg_mod = fes

            # separate observations by depth
            try:
                ind_deep = tg_mod.dataset.bathymetry.values > 200
            except:
                ind_deep = nemo1.dataset.bathymetry.values > 200
                print(f"Issue with bathy in {EXP}. Use bathy from GS1p1")

            try:
                X1, Y1 = obs.dataset[constit + "x"], tg_mod.dataset[constit + "x"]
                X2, Y2 = obs.dataset[constit + "y"], tg_mod.dataset[constit + "y"]

                plot_scatter_score(X1, Y1, X2, Y2, ind_deep, constit_str=constit,
                                   subtitle_str=[constit + "x (m)", constit + "y (m)"], title_str=EXP + "_complex")
            except:
                plt.close('all')
                plt.figure()
                plt.plot(0, 0);
                plt.title(EXP + ":" + constit)
                plt.savefig(config.dn_out + "PROCESSED/FIGS/" + constit + "_fail_" + EXP + "_complex.png")
                print(f"Issue with {constit} in {EXP}")

def plot_all_amp_pha_errors():
    # Loop over model runs and harmonic species.
    # Plot amp & pha with best fit stats against obs
    for EXP in ["ZPS_TIDE", "MES_TIDE", "FES2014"]:
        for constit in constit_list:
            if EXP == "ZPS_TIDE": tg_mod = nemo1
            if EXP == "MES_TIDE": tg_mod = nemo2
            if EXP == "FES2014": tg_mod = fes
            print(f"Plot {constit} amplitude and phase errors (deg): {EXP}")

            # separate observations by depth
            try:
                ind_deep = tg_mod.dataset.bathymetry.values > 200
            except:
                ind_deep = nemo1.dataset.bathymetry.values > 200
                print(f"Issue with bathy in {EXP}. Use bathy from GS1p1")

            try:
                tg_mod.dataset['A'], tg_mod.dataset['G'] = amp_pha_from_re_im(tg_mod.dataset[constit + "x"],
                                                                              tg_mod.dataset[constit + "y"])
                obs.dataset['A'], obs.dataset['G'] = amp_pha_from_re_im(obs.dataset[constit + "x"],
                                                                        obs.dataset[constit + "y"])

                X1, Y1 = obs.dataset.A, tg_mod.dataset.A
                X2, Y2 = modulo_align_phases_to_x(obs.dataset.G.values, tg_mod.dataset.G.values,
                                                  deg=True)  # preprocess w/ modulo arithmetic

                plot_scatter_score(X1, Y1, X2, Y2, ind_deep, constit_str=constit,
                                   subtitle_str=[constit + " Amplitude (m)", constit + " Phase (deg)"],
                                   title_str=EXP + "_amp_pha")

            except:
                plt.close('all')
                plt.figure()
                plt.plot(0, 0);
                plt.title(EXP + ":" + constit)
                plt.savefig(config.dn_out + "PROCESSED/FIGS/" + constit + "_fail_" + EXP + "_amp_pha.png")
                print(f"Issue with {constit} in {EXP}")

def plot_all_taylor_tides():
    # Loop over harmonic species.
    # Plot Taylor Tide diag of model and obs for each harmonic
    for constit in constit_list:

        count = 0
        for subset in ['shal', 'deep']:
            try:
                del II, z1obs, z2obs, z1mod, z2mod
            except:
                pass
            if subset == 'deep':
                # separate observations by depth
                II = nemo1.dataset.bathymetry.values > 200
            elif subset == 'shal':
                II = nemo1.dataset.bathymetry.values <= 200
            else:
                print(f"Not expecting that {subset}")

            z1obs, z2obs = obs.dataset[constit + 'x'][II], obs.dataset[constit + 'y'][II]
            if (0):
                # %%  Plot distributions of depth at observation locations
                plt.close('all')
                plt.figure()
                plt.plot(np.sort(obs.dataset.A[I].values))
                plt.xlabel('count')
                plt.ylabel(constit + ' Amp (m))')
                plt.title("distribution of A at observation sites")
                plt.legend()
                plt.savefig(config.dn_out + "PROCESSED/FIGS/dist_obsA_" + subset + ".png")

            # Obs
            if count == 0: # first pass. Obs: R=costheta=1, err=0
                R = np.array([1])  # R for obs
                rms_amp = np.array([rms_abs_amp(z1obs, z2obs)])
                rms_err = np.array([0])  # err for obs
            else: # second pass - other depth range. Obs: R=costheta=1, err=0
                R = np.hstack((R, 1))
                rms_amp = np.hstack((rms_amp, rms_abs_amp(z1obs, z2obs)))
                rms_err = np.hstack((rms_err, 0))

            # FES
            z1mod, z2mod = fes.dataset[constit + 'x'][II], fes.dataset[constit + 'y'][II]
            # obs_new, mod_new = tganalysis.match_missing_values(obs.dataset.M2x, tg_fes.dataset.M2x)
            # z1obs_new = obs_new.dataset.M2x
            # z1mod_new = mod_new.dataset.M2x
            # obs_new, mod_new = tganalysis.match_missing_values(obs.dataset.M2y, tg_fes.dataset.M2y)
            # z2obs_new = obs_new.dataset.M2y
            # z2mod_new = mod_new.dataset.M2y

            R = np.hstack((R, pearson_correl_coef(z1obs, z2obs, z1mod, z2mod)))
            rms_amp = np.hstack((rms_amp, rms_abs_amp(z1mod, z2mod)))
            rms_err = np.hstack((rms_err, rms_abs_error(z1obs, z2obs, z1mod, z2mod)))

            # nemo1
            del z1mod, z2mod
            z1mod, z2mod = nemo1.dataset[constit + 'x'][II], nemo1.dataset[constit + 'y'][II]
            # z1obs_new, z1mod_new = tganalysis.match_missing_values(obs.dataset.M2x, tg_mx2.dataset.M2x)
            # z2obs_new, z2mod_new = tganalysis.match_missing_values(obs.dataset.M2y, tg_mx2.dataset.M2y)

            R = np.hstack((R, pearson_correl_coef(z1obs, z2obs, z1mod, z2mod)))
            rms_amp = np.hstack((rms_amp, rms_abs_amp(z1mod, z2mod)))
            rms_err = np.hstack((rms_err, rms_abs_error(z1obs, z2obs, z1mod, z2mod)))

            # nemo2
            del z1mod, z2mod
            try:
                z1mod, z2mod = nemo2.dataset[constit + 'x'][II], nemo2.dataset[constit + 'y'][II]
            except:
                z1mod = np.nan
                z2mod = np.nan
            # z1obs_new, z1mod_new = tganalysis.match_missing_values(obs.dataset.M2x, tg_mv2.dataset.M2x)
            # z2obs_new, z2mod_new = tganalysis.match_missing_values(obs.dataset.M2y, tg_mv2.dataset.M2y)
            R = np.hstack((R, pearson_correl_coef(z1obs, z2obs, z1mod, z2mod)))
            rms_amp = np.hstack((rms_amp, rms_abs_amp(z1mod, z2mod)))
            rms_err = np.hstack((rms_err, rms_abs_error(z1obs, z2obs, z1mod, z2mod)))

            # TDISS-TEST
            del z1mod, z2mod
            try:
                z1mod, z2mod = tdiss.dataset[constit + 'x'][II], tdiss.dataset[constit + 'y'][II]
            except:
                z1mod = np.nan
                z2mod = np.nan
            # z1obs_new, z1mod_new = tganalysis.match_missing_values(obs.dataset.M2x, tg_mv2.dataset.M2x)
            # z2obs_new, z2mod_new = tganalysis.match_missing_values(obs.dataset.M2y, tg_mv2.dataset.M2y)
            R = np.hstack((R, pearson_correl_coef(z1obs, z2obs, z1mod, z2mod)))
            rms_amp = np.hstack((rms_amp, rms_abs_amp(z1mod, z2mod)))
            rms_err = np.hstack((rms_err, rms_abs_error(z1obs, z2obs, z1mod, z2mod)))

            count = count + 1

        label = ['obs:s', 'fes:s', 'gs1p1:s', 'gs1p7:s', 'tdiss:s',
                 'obs:d', 'fes:d', 'gs1p1:d', 'gs1p7:d', 'tdiss:d']
        npts = len(label)

        print(f"R= {[format(R[i], '.2f') for i in range(len(R))]}")
        print(f"rms_amp= {[format(rms_amp[i], '.2f') for i in range(len(R))]}")
        print(f"rms_err= {[format(rms_err[i], '.2f') for i in range(len(R))]}")
        print(label)

        ## Check cosine rule consistency
        A = rms_amp
        C = rms_err
        costheta = R

        B = rms_amp[0]
        for i in range(0, int(npts/2)):
            print(
                f"{constit}:{label[i]}: sqrt(A^2+B^2-2ABcos(theta))={np.sqrt(A[i] ** 2 + B ** 2 - 2 * A[i] * B * costheta[i])}. C={C[i]}")
        B = rms_amp[int(npts/2)]
        for i in range(int(npts/2), int(npts)):
            print(
                f"{constit}:{label[i]}: sqrt(A^2+B^2-2ABcos(theta))={np.sqrt(A[i] ** 2 + B ** 2 - 2 * A[i] * B * costheta[i])}. C={C[i]}")

        # Create TaylorTide plot template
        tt = TaylorTide(
            r_obs=rms_amp[0],
            rms_amp_max=0.7,
            rms_amp_contours=[0.2, 0.4, 0.6],
            rms_err_contours=[0.2, 0.4, 0.6],
            cos_theta_lines=[0.3, 0.6, 0.9],
        )
        # Add data to axes
        tt.ax.scatter(rms_amp[1:int(npts/2)] * R[1:int(npts/2)],
                      rms_amp[1:int(npts/2)] * np.sqrt(1 - R[1:int(npts/2)] ** 2),
                      s=10, c=['r', 'k', 'g', 'y'])
        # manual legend
        colors = ['red', 'black', 'green', 'yellow']
        lines = [Line2D([0], [0], color=c, linewidth=3, linestyle='dotted') for c in colors]
        labels = ["FES2014", "GS1p1", "GS1p7", "TDISS-TEST"]
        plt.legend(lines, labels)

        plt.title(constit + ':shallow' + " N=" + str(int(np.isfinite(z1obs).sum())))
        plt.savefig(config.dn_out + "PROCESSED/FIGS/Taylor_" + constit + "_shallow.png")

        # Create TaylorTide plot template
        tt = TaylorTide(
            r_obs=rms_amp[4],
            rms_amp_max=0.61,
            rms_amp_contours=[0.2, 0.4, 0.6],
            rms_err_contours=[0.2, 0.4, 0.6],
            cos_theta_lines=[0.3, 0.6, 0.9],
        )
        # Add data to axes
        tt.ax.scatter(rms_amp[int(npts/2)+1::] * R[int(npts/2)+1::],
                      rms_amp[int(npts/2)+1::] * np.sqrt(1 - R[int(npts/2)+1::] ** 2),
                      s=10, c=['r', 'k', 'g', 'y'])
        # manual legend
        colors = ['red', 'black', 'green', 'yellow']
        lines = [Line2D([0], [0], color=c, linewidth=3, linestyle='dotted') for c in colors]
        labels = ["FES2014", "GS1p1", "GS1p7", "TDISS-TEST"]
        plt.legend(lines, labels)

        plt.title(constit + ':deep' + " N=" + str(int(np.isfinite(z1obs).sum())))
        plt.savefig(config.dn_out + "PROCESSED/FIGS/Taylor_" + constit + "_deep.png")

def plot_overlay_taylor_tides():
    # Loop over harmonic species.
    # Plot Taylor Tide diag of model and obs for each harmonic. Overlay on two (deep/shallow) plots as trees

    nsim=4 # number of simulations + obs.  labels = ["GS1p1", "GS1p2", "FES2014", "obs"]

    #for subset in ['shal']:
    for subset in ["shal", "deep"]:
      #for constit_family_list in [["M2", "S2", "N2", "K2"]]:
      for constit_family_list in [["M2", "S2", "N2", "K2"], ["O1", "Q1", "P1"]]:
        if "M2" in constit_family_list: family_str = "semi-diurnal"
        if "O1" in constit_family_list: family_str = "diurnal"
        R = np.zeros((nsim, len(constit_family_list)))
        rms_amp = np.zeros((nsim, len(constit_family_list)))
        rms_err = np.zeros((nsim, len(constit_family_list)))
        label = {}
        for count, constit in enumerate(constit_family_list):
            try:
                del II, z1obs, z2obs, z1mod, z2mod
            except:
                pass
            if subset == "deep":
                # separate observations by depth
                II = nemo1.dataset.bathymetry.values > 200
            elif subset == "shal":
                II = nemo1.dataset.bathymetry.values <= 200
            else:
                print(f"Not expecting that {subset}")

            z1obs, z2obs = obs.dataset[constit + 'x'][II], obs.dataset[constit + 'y'][II]

            # Obs
            #if count == 0:
            #    R = np.array([1])  # R for obs
            #    rms_amp = np.array([rms_abs_amp(z1obs, z2obs)])
            #    rms_err = np.array([0])  # err for obs
            #    label = [constit]
            #else:
            #    R = np.hstack((R, 1))
            #    rms_amp = np.hstack((rms_amp, rms_abs_amp(z1obs, z2obs)))
            #    rms_err = np.hstack((rms_err, 0))
            #    label.append(constit)

            R[0,count] = 1
            rms_amp[0,count] = rms_abs_amp(z1obs, z2obs)
            rms_err[0,count] = 0
            label[0,count] = 'obs:'+constit


            # nemo1
            try:
                del z1mod, z2mod
            except:
                pass
            z1mod, z2mod = nemo1.dataset[constit + 'x'][II], nemo1.dataset[constit + 'y'][II]

            R[1,count] = pearson_correl_coef(z1obs, z2obs, z1mod, z2mod)
            rms_amp[1,count] = rms_abs_amp(z1mod, z2mod)
            rms_err[1,count] = rms_abs_error(z1obs, z2obs, z1mod, z2mod)
            label[1,count] = 'ZPS_TIDE:'+constit



            # nemo2
            del z1mod, z2mod
            #try:
            z1mod, z2mod = nemo2.dataset[constit + 'x'][II], nemo2.dataset[constit + 'y'][II]
            #except:
            #    z1mod = np.nan
            #    z2mod = np.nan

            R[2,count] = pearson_correl_coef(z1obs, z2obs, z1mod, z2mod)
            rms_amp[2,count] = rms_abs_amp(z1mod, z2mod)
            rms_err[2,count] = rms_abs_error(z1obs, z2obs, z1mod, z2mod)
            label[2,count] = 'MES_TIDE:'+constit


            # FES
            del z1mod, z2mod
            z1mod, z2mod = fes.dataset[constit + 'x'][II], fes.dataset[constit + 'y'][II]

            R[3,count] = pearson_correl_coef(z1obs, z2obs, z1mod, z2mod)
            rms_amp[3,count] = rms_abs_amp(z1mod, z2mod)
            rms_err[3,count] = rms_abs_error(z1obs, z2obs, z1mod, z2mod)
            label[3,count] = 'FES:'+constit


        print(subset)
        #print(f"R= {[format(R[i], '.2f') for i in range(len(R))]}")
        #print(f"rms_amp= {[format(rms_amp[i], '.2f') for i in range(len(R))]}")
        #print(f"rms_err= {[format(rms_err[i], '.2f') for i in range(len(R))]}")
        print(label)

        ## Check cosine rule consistency. Output data.
        for count, constit in enumerate(constit_family_list):
        #for i in range(len(constit_family_list)):
            B = rms_amp[0,count]
            A = rms_amp[1:nsim,count]
            C = rms_err[1:nsim,count]
            costheta = R[1:nsim,count]

            print("Check cosine rule consistency")
            for j in range(0,nsim-1): # model runs: fes, nemo1, gs1p2, tdiss
                print(
                    f"{label[1+j,count]}: sqrt(A^2+B^2-2ABcos(theta))={np.sqrt(A[j] ** 2 + B ** 2 - 2 * A[j] * B * costheta[j])}. C={C[j]}")

            print("Output table of data")
            for j in range(0, nsim - 1):  # model runs: fes, nemo1, gs1p2, tdiss
                print(
                    f"{label[1 + j, count]:11s}: (A, B, C, theta)={A[j]:.3f}, {B:.3f}, {C[j]:.3f}, {np.arccos(costheta[j]) * 180 / np.pi:.1f}")
            del B, A, C, costheta



        # Create TaylorTide plot template
        if family_str == "semi-diurnal":
            fig_amp_max = 0.6
            if subset == "shal":   fig_err_contours = [0.1, 0.2, 0.3]
            elif subset == "deep": fig_err_contours = [0.05, 0.1, 0.15]

        elif family_str == "diurnal":
            fig_amp_max = 0.1
            if subset == "shal":   fig_err_contours = [0.01, 0.02, 0.03]
            elif subset == "deep": fig_err_contours = [0.01, 0.02, 0.03]



        tt = TaylorTide(
            r_obs=rms_amp[0,0],
            rms_amp_max=fig_amp_max,
            rms_amp_contours=[],
            #rms_amp_contours=[0.2, 0.4, 0.6],
            rms_err_contours=fig_err_contours,
            cos_theta_lines=[0.3, 0.6, 0.9],
            err_contour_flag=True,
        )

        ## Loop over constituents to add data
        for i in range(len(constit_family_list)):
            # Add data to axes
            tt.ax.scatter(rms_amp[1:nsim,i] * R[1:nsim,i],
                          rms_amp[1:nsim,i] * np.sqrt(1 - R[1:nsim,i] ** 2),
                          s=60, c=[c_ZPS_TIDE, c_MES_TIDE, c_FES], zorder=10, clip_on=False)
                            #s = 20, c = ['b', 'g', 'r', 'y'], zorder = 10, clip_on = False)
            # Add vectors between points and obs
            tt.ax.plot([np.repeat(rms_amp[0,i],nsim-1), rms_amp[1:nsim,i] * R[1:nsim,i]],
                       [np.zeros(nsim-1), rms_amp[1:nsim,i] * np.sqrt(1 - R[1:nsim,i] ** 2)])
            if nsim != 4: print('Colours and lines not as expected here')
            #tt.ax.lines[-4].set_color('k')
            tt.ax.lines[-3].set_color(c_ZPS_TIDE) #'g')
            tt.ax.lines[-2].set_color(c_MES_TIDE) #'r')
            tt.ax.lines[-1].set_color(c_FES) #'y')

            # Place constituent labels
            #tt.ax.text( rms_amp[0,i], -0.036*fig_amp_max, constit_family_list[i], rotation=0, color='b')
            xpos = rms_amp[0,i]
            ypos = 0.036*fig_amp_max
            if constit_family_list[i]=="M2": xpos = rms_amp[0,i] - 0.07
            if constit_family_list[i]=="K2": xpos = rms_amp[0,i] - 0.03; ypos = ypos + 0.02
            if constit_family_list[i]=="N2": xpos = rms_amp[0,i] - 0.02; ypos = ypos + 0.01
            if constit_family_list[i]=="Q1": xpos = rms_amp[0,i] - 0.012
            if constit_family_list[i]=="O1": xpos = rms_amp[0,i] + 0.002
            if constit_family_list[i]=="P1": xpos = rms_amp[0,i] + 0.002
            tt.ax.text( xpos, ypos, constit_family_list[i], fontsize=12,  rotation=0, color='b')

            # Draw obs dot (smaller)
            tt.ax.scatter(rms_amp[0,i], 0, s=40, color=c_obs, zorder=20, clip_on=False)  # obs point

        # Add indicative timing labels
        if family_str == "diurnal":
            tt.ax.text( fig_amp_max*np.cos(np.pi/6), fig_amp_max*np.sin(np.pi/6), "      (2hr)", color='k')
            tt.ax.text( fig_amp_max*np.cos(np.pi/24), fig_amp_max*np.sin(np.pi/24), " (30min)", color='k')
        elif family_str == "semi-diurnal":
            tt.ax.text( fig_amp_max*np.cos(np.pi/6), fig_amp_max*np.sin(np.pi/6),  "      (1hr)", color='k')
            tt.ax.text( fig_amp_max*np.cos(np.pi/24), fig_amp_max*np.sin(np.pi/24)," (15min)", color='k')

        # manual legend (once)
        if subset == "shal" and family_str == "semi-diurnal":
            #colors = [ 'green', 'red', 'yellow', 'blue']
            colors = [ c_ZPS_TIDE, c_MES_TIDE, c_FES, c_obs]
            sizes = [6, 6, 6, 4]
            lines = [Line2D([], [], color=colors[ci], markersize=sizes[ci], marker='o', linestyle='None') for ci in range(len(colors))]
            labels = ["ZPS_TIDE", "MES_TIDE", "FES2014", "observations"]
            plt.legend(lines, labels, loc='upper left', framealpha=1)

        #plt.title(subset + ":" + family_str)
        if subset == "shal": sym = "<"
        elif subset == "deep": sym = ">"
        plt.title(family_str + ": water depth " + sym + " 200m")
        plt.savefig(config.dn_out + "PROCESSED/FIGS/Taylor_" + "_" + subset + "_" + family_str + "_tree.png")


        # Plot Taylor with polars
        #########################
        if(0):
            #theta = np.linspace(0, 1, 100)* 45
            #r = np.ones((100,1))
            r = np.arange(0, 1, 0.01)
            theta = 0.25 * np.pi * r

            fig, ax_polar = plt.subplots(subplot_kw={'projection': 'polar'})
            ax_polar.grid(False)

            ## Loop over constituents to create plot
            for i in range(len(constit_family_list)):
                # Add data to axes
                ax_polar.scatter( np.arccos(R[0:nsim,i]), rms_amp[0:nsim,i],
                              s=20, c=['b', 'g', 'r', 'y'])
                # Add vectors between points and obs
                ax_polar.plot([np.zeros(nsim-1), np.arccos(R[1:nsim,i])],
                              [np.repeat(rms_amp[0,i],nsim-1), rms_amp[1:nsim,i]])
                #ax_polar.plot([np.repeat(rms_amp[0,i],nsim-1), rms_amp[1:nsim,i] * R[1:nsim,i]],
                #           [np.zeros(nsim-1), rms_amp[1:nsim,i] * np.sqrt(1 - R[1:nsim,i] ** 2)])
                if nsim != 4: print('Colours and lines not as expected here')
                #tt.ax.lines[-4].set_color('k')
                ax_polar.lines[-3].set_color('g')
                ax_polar.lines[-2].set_color('r')
                ax_polar.lines[-1].set_color('y')

                # Blue arc through obs and label
                ax_polar.plot(theta,
                                 np.repeat(rms_amp[0,i],100), '-', color='blue')  # obs point
                ax_polar.text( 0, rms_amp[0,i]+0.02, constit_family_list[i], rotation=0, color='b')

                # Add error arc for each constituent
                #theta_rad = np.array([0, np.pi/4, np.pi/2, np.pi*0.75, np.pi]) #
                theta_rad = np.arange(0, np.pi, 0.01)
                for ir in np.linspace(rms_err[1,i]/5, rms_err[1,i], 4):
                    rr = ir # rms_err[1,i]
                    xx = rms_amp[0,i] + rr * np.cos(theta_rad)
                    yy = rr * np.sin(theta_rad)
                    rrr = np.sqrt( xx**2 + yy**2 )
                    phi = np.arctan2(yy,xx)
                    ax_polar.plot( phi, rrr, '--', color='grey')
            ax_polar.set_thetamin(0)
            ax_polar.set_thetamax(35)
            ax_polar.set_rmin(0.45)
            ax_polar.set_rmax(0.6)
            ax_polar.set_rorigin(0)

            #plt.show()
            plt.savefig(config.dn_out + "PROCESSED/FIGS/Taylor_" + "_" + subset + "_" + family_str + "_tree_polar.png")

def make_table():
    nsim=4 # number of simulations + obs.  labels = ["GS1p1", "GS1p2", "FES2014", "obs"]
    dict = {'shal':{"M2":{"obs":{}, "ZPS_TIDE":{}, "MES_TIDE":{}, "FES":{}},
                    "S2":{"obs":{}, "ZPS_TIDE":{}, "MES_TIDE":{}, "FES":{}},
                    "N2":{"obs":{}, "ZPS_TIDE":{}, "MES_TIDE":{}, "FES":{}},
                    "K2":{"obs":{}, "ZPS_TIDE":{}, "MES_TIDE":{}, "FES":{}},
                    "O1":{"obs":{}, "ZPS_TIDE":{}, "MES_TIDE":{}, "FES":{}},
                    "Q1":{"obs":{}, "ZPS_TIDE":{}, "MES_TIDE":{}, "FES":{}},
                    "P1":{"obs":{}, "ZPS_TIDE":{}, "MES_TIDE":{}, "FES":{}}},
            'deep':{"M2":{"obs":{}, "ZPS_TIDE":{}, "MES_TIDE":{}, "FES":{}},
                    "S2":{"obs":{}, "ZPS_TIDE":{}, "MES_TIDE":{}, "FES":{}},
                    "N2":{"obs":{}, "ZPS_TIDE":{}, "MES_TIDE":{}, "FES":{}},
                    "K2":{"obs":{}, "ZPS_TIDE":{}, "MES_TIDE":{}, "FES":{}},
                    "O1":{"obs":{}, "ZPS_TIDE":{}, "MES_TIDE":{}, "FES":{}},
                    "Q1":{"obs":{}, "ZPS_TIDE":{}, "MES_TIDE":{}, "FES":{}},
                    "P1":{"obs":{}, "ZPS_TIDE":{}, "MES_TIDE":{}, "FES":{}}}}
    for subset in ["shal", "deep"]:
      for constit_family_list in [["M2", "S2", "N2", "K2"], ["O1", "Q1", "P1"]]:

        if "M2" in constit_family_list: family_str = "semi-diurnal"
        if "O1" in constit_family_list: family_str = "diurnal"
        R = np.zeros((nsim, len(constit_family_list)))
        rms_amp = np.zeros((nsim, len(constit_family_list)))
        rms_err = np.zeros((nsim, len(constit_family_list)))
        label = {}
        for count, constit in enumerate(constit_family_list):
            try:
                del II, z1obs, z2obs, z1mod, z2mod
            except:
                pass
            if subset == "deep":
                # separate observations by depth
                II = nemo1.dataset.bathymetry.values > 200
            elif subset == "shal":
                II = nemo1.dataset.bathymetry.values <= 200
            else:
                print(f"Not expecting that {subset}")

            z1obs, z2obs = obs.dataset[constit + 'x'][II], obs.dataset[constit + 'y'][II]

            # Obs
            R[0,count] = 1
            rms_amp[0,count] = rms_abs_amp(z1obs, z2obs)
            rms_err[0,count] = 0
            label[0,count] = 'obs:'+constit

            # Store information
            dict[subset][constit]['obs'] = {'amp': f"{rms_amp[0,count]:.3f}",
                                              'theta': 0,
                                              'error': 0}


            # nemo1
            try:
                del z1mod, z2mod
            except:
                pass
            z1mod, z2mod = nemo1.dataset[constit + 'x'][II], nemo1.dataset[constit + 'y'][II]

            R[1,count] = pearson_correl_coef(z1obs, z2obs, z1mod, z2mod)
            rms_amp[1,count] = rms_abs_amp(z1mod, z2mod)
            rms_err[1,count] = rms_abs_error(z1obs, z2obs, z1mod, z2mod)
            label[1,count] = 'ZPS_TIDE:'+constit

            # Store information
            dict[subset][constit]['ZPS_TIDE'] = {'amp': f"{rms_amp[1,count]:.3f}",
                                                 'theta': f"{np.arccos(R[1,count])*180/np.pi:.1f}",
                                                 'error': f"{rms_err[1,count]:.3f}"}


            # nemo2
            del z1mod, z2mod
            #try:
            z1mod, z2mod = nemo2.dataset[constit + 'x'][II], nemo2.dataset[constit + 'y'][II]
            #except:
            #    z1mod = np.nan
            #    z2mod = np.nan

            R[2,count] = pearson_correl_coef(z1obs, z2obs, z1mod, z2mod)
            rms_amp[2,count] = rms_abs_amp(z1mod, z2mod)
            rms_err[2,count] = rms_abs_error(z1obs, z2obs, z1mod, z2mod)
            label[2,count] = 'MES_TIDE:'+constit

            # Store information
            dict[subset][constit]['MES_TIDE'] = {'amp': f"{rms_amp[2,count]:.3f}",
                                                 'theta': f"{np.arccos(R[2,count])*180/np.pi:.1f}",
                                                 'error': f"{rms_err[2,count]:.3f}"}


            # FES
            del z1mod, z2mod
            z1mod, z2mod = fes.dataset[constit + 'x'][II], fes.dataset[constit + 'y'][II]

            R[3,count] = pearson_correl_coef(z1obs, z2obs, z1mod, z2mod)
            rms_amp[3,count] = rms_abs_amp(z1mod, z2mod)
            rms_err[3,count] = rms_abs_error(z1obs, z2obs, z1mod, z2mod)
            label[3,count] = 'FES:'+constit

            # Store information
            dict[subset][constit]['FES'] = {'amp': f"{rms_amp[3,count]:.3f}",
                                                 'theta': f"{np.arccos(R[3,count])*180/np.pi:.1f}",
                                                 'error': f"{rms_err[3,count]:.5f}"}
    print(dict)

    def dict_to_latex(data):
        # lay out the constituents as rows. Cell: "amp (err, phase)"
        latex_table = "\\begin{tabular}{|l|c|c|c|c|}\n"
        latex_table += "\\hline\n"
        latex_table += " & obs & ZPS\\_TIDE & MES\\_TIDE & FES \\\\\n"
        latex_table += "\\hline\n"
        for depth_category, fields in data.items():
            for tide_con, values in fields.items():
                latex_table += f"{tide_con} & "
                latex_table += " & ".join(
                    [f"{values[method]['amp']} ({values[method]['error']}, {values[method]['theta']})" for method in values])
                latex_table += " \\\\\n"
            latex_table += "\\hline\n"
        latex_table += "\\end{tabular}"

        # lay out the sea level rms by dataset (column). Sum of errors
        latex_table += "\n"
        latex_table += "\\begin{tabular}{|l|c|c|c|}\n"
        latex_table += "\\hline\n"
        latex_table += " & ZPS\\_TIDE & MES\\_TIDE & FES \\\\\n"
        latex_table += "\\hline\n"
        for depth_category in ["shal", "deep"]:

            val = np.array([float(dict[depth_category][c]['ZPS_TIDE']['error']) for c in ['M2', 'S2', 'N2', 'K2', 'O1', 'Q1', 'P1']]).sum()
            latex_table += f"$\\eta$ & {val/np.sqrt(2):.3f} & "
            val = np.array([float(dict[depth_category][c]['MES_TIDE']['error']) for c in ['M2', 'S2', 'N2', 'K2', 'O1', 'Q1', 'P1']]).sum()
            latex_table += f"{val/np.sqrt(2):.3f} & "
            val = np.array([float(dict[depth_category][c]['FES']['error']) for c in ['M2', 'S2', 'N2', 'K2', 'O1', 'Q1', 'P1']]).sum()
            latex_table += f"{val/np.sqrt(2):.3f} \\\\"
            latex_table += "\\hline\n"

        latex_table += "\\end{tabular}"

        return latex_table

    latex_table = dict_to_latex(dict)
    print(latex_table)
    with open('output.tex', 'w') as file:
        file.write(latex_table)

def plot_cloud():
    # Taylor Tide with cloud of all data points
    """
    This is not quite right. Or maybe I need to document it better.
    The coloured dots should be mean model error and mean obs _error_. The black
     dots should be model err.
     Not sure why there is a magenta and a red dot.
     The count variable doesn't seem to be doing the right thing by not changing
     """
    #for constit in constit_list:
    constit = "M2"
    EXP =  "ZPS_TIDE" #, "MES_TIDE", "FES2014",
    if EXP == "ZPS_TIDE": tg_mod = nemo1
    if EXP == "MES_TIDE": tg_mod = nemo2
    if EXP == "FES2014": tg_mod = fes


    if(1):
        count = 0
        subset = "shal"
        try:
            del II, z1obs, z2obs, z1mod, z2mod, rms_amp
        except:
            pass
        if subset == 'deep':
            # separate observations by depth
            II = nemo1.dataset.bathymetry.values > 200
        elif subset == 'shal':
            II = nemo1.dataset.bathymetry.values <= 200
        else:
            print(f"Not expecting that {subset}")

        z1obs, z2obs = obs.dataset[constit + 'x'][II], obs.dataset[constit + 'y'][II]

        # Obs
        if count == 0:
            R = np.array([1])  # R for obs
            rms_amp     = np.array([rms_abs_amp(z1obs, z2obs)])
            rms_amp_obs = np.array([rms_abs_amp(z1obs, z2obs)])
            rms_err = np.array([0])  # err for obs


        z1mod, z2mod = tg_mod.dataset[constit + 'x'][II], tg_mod.dataset[constit + 'y'][II]

        # model - averaged
        if count == 1:
            R = np.hstack((R, pearson_correl_coef(z1obs, z2obs, z1mod, z2mod)))
            rms_amp     = np.hstack((rms_amp,     rms_abs_amp(z1mod, z2mod)))
            rms_amp_obs = np.hstack((rms_amp_obs, rms_abs_amp(z1obs, z2obs)))
            rms_err = np.hstack((rms_err, rms_abs_error(z1obs, z2obs, z1mod, z2mod)))

        n_length = len(z1obs)
        for ii in range(n_length):
            R = np.hstack((R, pearson_correl_coef(z1obs[ii], z2obs[ii], z1mod[ii], z2mod[ii])))
            rms_amp     = np.hstack((rms_amp,     rms_abs_amp(z1mod[ii], z2mod[ii])))
            rms_amp_obs = np.hstack((rms_amp_obs, rms_abs_amp(z1obs[ii], z2obs[ii])))
            rms_err = np.hstack((rms_err, rms_abs_error(z1obs[ii], z2obs[ii], z1mod[ii], z2mod[ii])))




        #label = ['obs:s', 'fes:s', 'gs1p1:s', 'gs1p2:s',
        #         'obs:d', 'fes:d', 'gs1p1:d', 'gs1p2:d']

        print(f"R= {[format(R[i], '.2f') for i in range(len(R))]}")
        print(f"rms_amp= {[format(rms_amp[i], '.2f') for i in range(len(R))]}")
        print(f"rms_err= {[format(rms_err[i], '.2f') for i in range(len(R))]}")
        #print(label)

        ## Check cosine rule consistency
        A = rms_amp
        C = rms_err
        costheta = R

        B = rms_amp_obs #rms_amp[0]
        for i in range(0, 4):
            print(
                f"{i}: sqrt(A^2+B^2-2ABcos(theta))={np.sqrt(A[i] ** 2 + B[i] ** 2 - 2 * A[i] * B[i] * costheta[i])}. C={C[i]}")
        plt.figure()
        plt.plot( np.sqrt(A ** 2 + B ** 2 - 2 * A * B * costheta), C, '+')
        plt.xlabel("sqrt(A^2+B^2-2ABcos(theta))")
        plt.ylabel("C")
        plt.title(f"{EXP}: Check cosine rule consistency for {constit}")
        plt.show()

        # Create TaylorTide plot template
        tt = TaylorTide(
            r_obs=rms_amp[0],
            rms_amp_max=2, #1.75,
            rms_amp_contours=[0.2, 0.4, 0.6, 1, 1.5],
            rms_err_contours=[0.2, 0.4, 0.6],
            cos_theta_lines=[0.3, 0.6, 0.9],
        )
        # Add data to axes
        tt.ax.scatter(rms_amp[2:n_length] * R[2:n_length], rms_amp[2:n_length] * np.sqrt(1 - R[2:n_length] ** 2), s=10, c='k')
        tt.ax.scatter(rms_amp_obs[2:n_length], -np.ones_like(rms_amp_obs)[2:n_length]*0.02, s=10, c='b')
        tt.ax.scatter(rms_amp[0] * R[0], rms_amp[0] * np.sqrt(1 - R[0] ** 2), s=20, c='b')  # obs ave
        tt.ax.scatter(rms_amp[1] * R[1], rms_amp[1] * np.sqrt(1 - R[1] ** 2), s=20, c='r')  # model ave
        #tt.ax.scatter(rms_amp[2] * R[2], rms_amp[2] * np.sqrt(1 - R[2] ** 2), s=20, c='m')

        for i in [0,1,2]:
            print(f"{rms_amp[i] * R[i], rms_amp[i] * np.sqrt(1 - R[i]**2)}")

        plt.title(EXP + " " + constit + ':dots')
        plt.savefig(config.dn_out + "PROCESSED/FIGS/Taylor_" + constit + "_" + EXP + "_shallow_SPECIAL.png")

def plot_east_coast_usa():
    ylims = [35, 45]  # [-80,80] #[45,60]
    xlims = [-76, -60]  # [-180,180] #[-60, -25]

    plt.close('all')
    #plt.figure()
    fig, axs = plt.subplots(3, 2, figsize=(10, 5), subplot_kw={'projection': ccrs.PlateCarree()})
    axs[0,0].coastlines()
    axs[0,0].scatter(obs.dataset.longitude, obs.dataset.latitude, c=obs.dataset['A'], s=20,
                     vmin=0, vmax=0.45, cmap='Spectral')
    axs[0,0].set_title('obs')

    axs[0,1].coastlines()
    axs[0,1].scatter(fes.dataset.longitude, fes.dataset.latitude, c=fes.dataset['A'], s=20,
                     vmin=0, vmax=0.45, cmap='Spectral')
    axs[0,1].set_title('FES')


    axs[1,0].coastlines()
    axs[1,0].scatter(nemo1.dataset.longitude, nemo1.dataset.latitude, c=nemo1.dataset['A'], s=20,
                     vmin=0, vmax=0.45, cmap='Spectral')
    axs[1,0].set_title('GS1P1')


    axs[1,1].coastlines()
    im = axs[1,1].scatter(gs1p2.dataset.longitude, gs1p2.dataset.latitude, c=gs1p2.dataset['A'], s=20,
                     vmin=0, vmax=0.45, cmap='Spectral')
    axs[1,1].set_title('GS1P2')

    plt.setp(axs, xlim=xlims, ylim=ylims)

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)

    plt.suptitle("M2 amp")
    # plt.show()
    plt.savefig(config.dn_out + "PROCESSED/FIGS/scatter_amp_on_map.png")

def load_obs_data():
    # Load data as tidegauge objects
    ## Harmonise definitions: negate M2y (and phase) in NEMO - done in preprocessing.
    ################################
    obs = coast.Tidegauge(dataset=xr.open_dataset(config.dn_out+"PROCESSED/obs_extracted.nc"))
    obs.dataset['A'], obs.dataset['G'] = amp_pha_from_re_im(obs.dataset.M2x, obs.dataset.M2y)
    obs.dataset['longitude'] = ((obs.dataset.longitude%360) + 180)%360 - 180  # set longitude : [-180,180]
    obs.dataset['G'] = (obs.dataset.G+180)%360-180  # set phase: -180,180
    return obs

def load_fes_data():
    # Load FES data as tidegauge objects
    ## Harmonise definitions: negate M2y (and phase) in NEMO - done in preprocessing.
    ################################
    fes = coast.Tidegauge(dataset=xr.open_dataset(config.dn_out+"PROCESSED/FES2014_extracted.nc"))
    fes.dataset['A'], fes.dataset['G'] = amp_pha_from_re_im(fes.dataset.M2x, fes.dataset.M2y)
    fes.dataset['G'] = (fes.dataset.G+180)%360-180  # set phase: -180,180
    return fes        

def load_nemo_data(filename="GS1p7_triads_extracted.nc"):

    nemo = coast.Tidegauge(dataset=xr.open_dataset(config.dn_out+"PROCESSED/"+filename))
    nemo.dataset['A'], nemo.dataset['G'] = amp_pha_from_re_im(nemo.dataset.M2x, nemo.dataset.M2y)

    nemo.dataset['G'] = +(nemo.dataset.G+180)%360-180  # set phase: -180,180
    #nemo.dataset['M2y'] = +nemo.dataset.M2y
    nemo.dataset = nemo.dataset.drop_dims(["nvertex", "z_dim"])  # drop unwanted dimensions and associated variables
    return nemo



def load_data(filename1:str="GS1p1_tide_final_extracted.nc", filename2:str="GS1p7_triads_extracted.nc"):
    """ 
    Loads the observational and FES2014 data
    Also loads two additional filenames
    """
    # Load data as tidegauge objects
    ## Harmonise definitions: negate M2y (and phase) in NEMO - done in preprocessing.
    ################################

    obs = load_obs_data()
    fes = load_fes_data()

    nemo1 = load_nemo_data(filename1)
    nemo2 = load_nemo_data(filename2)


    #tdiss = coast.Tidegauge(dataset=xr.open_dataset(config.dn_out+"PROCESSED/G1p2_full_extracted.nc"))
    #tdiss.dataset['A'], tdiss.dataset['G'] = amp_pha_from_re_im(tdiss.dataset.M2x, tdiss.dataset.M2y)

    #tdiss.dataset['G'] = (tdiss.dataset.G+180)%360-180  # set phase: -180,180
    #tdiss.dataset['M2y'] = +tdiss.dataset.M2y
    #tdiss.dataset = tdiss.dataset.drop_dims(["nvertex", "z_dim"])  # drop unwanted dimensions and associated variables
    if(0):
        nemo1 = coast.Tidegauge(dataset=xr.open_dataset(config.dn_out+"PROCESSED/GS1p1_tide_final_extracted.nc"))
        nemo1.dataset['A'], nemo1.dataset['G'] = amp_pha_from_re_im(nemo1.dataset.M2x, nemo1.dataset.M2y)

        nemo1.dataset['G'] = (nemo1.dataset.G+180)%360-180  # set phase: -180,180
        nemo1.dataset['M2y'] = +nemo1.dataset.M2y
        nemo1.dataset = nemo1.dataset.drop_dims(["nvertex", "nav_lev", "z_dim"])  # drop unwanted dimensions and associated variables

        nemo2 = coast.Tidegauge(dataset=xr.open_dataset(config.dn_out+"PROCESSED/GS1p7_triads_extracted.nc"))
        #gs1p2.dataset = gs1p2.dataset.drop_vars(["G", "A"]) # Unlabelled amp/pha parts should not be in the file
        nemo2.dataset['A'], nemo2.dataset['G'] = amp_pha_from_re_im(nemo2.dataset.M2x, nemo2.dataset.M2y)

        nemo2.dataset['G'] = +(nemo2.dataset.G+180)%360-180  # set phase: -180,180
        nemo2.dataset['M2y'] = +nemo2.dataset.M2y
        nemo2.dataset = nemo2.dataset.drop_dims(["nvertex", "z_dim"])  # drop unwanted dimensions and associated variables

    return obs, nemo1, nemo2, fes

def scatter_phase_on_map(list_of_datasets, list_of_titles):
    """ Primary purpose: Check that the phases are all aligned """
    ## Basic map plot of phases
    # Phase alignment

    #extent = [90, 125, -10, 10]
    #central_lat = np.mean(extent[2:])
    #central_lon = np.mean(extent[:2])

    ylims = [-80,80] #[45,60]
    xlims = [-180,180] #[-60, -25]

    plt.close('all')

    fig, ax = plt.subplots(2, 2, figsize=(10,5),
                        subplot_kw={'projection': ccrs.PlateCarree()})
    #ax = plt.axes(projection=ccrs.AlbersEqualArea(central_lon, central_lat))
    #ax.set_extent(extent)

    ax[0,0].coastlines(resolution='10m', color='k', linestyle='-', alpha=1)
    tt = ax[0,0].scatter(list_of_datasets[0].dataset.longitude, list_of_datasets[0].dataset.latitude,
                         c=list_of_datasets[0].dataset['G'], s=10, cmap='Spectral')
    #plt.colorbar()
    ax[0,0].set_title(list_of_titles[0]) #('obs')
    ax[0,0].set_xlim(xlims)
    ax[0,0].set_ylim(ylims)

    ax[0,1].coastlines(resolution='10m', color='k', linestyle='-', alpha=1)
    ax[0,1].scatter(list_of_datasets[1].dataset.longitude, list_of_datasets[1].dataset.latitude,
                     c=list_of_datasets[1].dataset['G'], s=10, cmap='Spectral')
    #plt.colorbar()
    ax[0,1].set_title(list_of_titles[1])
    ax[0,1].set_xlim(xlims)
    ax[0,1].set_ylim(ylims)

    ax[1,0].coastlines(resolution='10m', color='k', linestyle='-', alpha=1)
    ax[1,0].scatter(list_of_datasets[2].dataset.longitude, list_of_datasets[2].dataset.latitude,
                     c=list_of_datasets[2].dataset['G'], s=10, cmap='Spectral')
    #plt.colorbar()
    ax[1,0].set_title(list_of_titles[2])
    ax[1,0].set_xlim(xlims)
    ax[1,0].set_ylim(ylims)

    ax[1,1].coastlines(resolution='10m', color='k', linestyle='-', alpha=1)
    ax[1,1].scatter(list_of_datasets[3].dataset.longitude, list_of_datasets[3].dataset.latitude,
                     c=list_of_datasets[3].dataset['G'], s=10, cmap='Spectral')
    #plt.colorbar()
    ax[1,1].set_title(list_of_titles[3])
    ax[1,1].set_xlim(xlims)
    ax[1,1].set_ylim(ylims)

    cbar_ax = fig.add_axes([0.1, 0.1, 0.8, 0.02])
    plt.colorbar(tt , orientation='horizontal', cax=cbar_ax)
    plt.suptitle("M2 phase")
    #plt.show()
    fig.savefig(config.dn_out+"PROCESSED/FIGS/scatter_M2_phase_on_map.png")
    plt.close('all')

### Main script
# load data
obs, nemo1, nemo2, fes = load_data(filename1="GS1p1_tide_final_extracted.nc", filename2="GS1p7_triads_extracted.nc")

# Compare maps
if(0):
 # Plot utility for COAsT.tidegauge objects
 obs.plot_on_map_multiple([obs], color_var_str="G")
 #plt.savefig(config.dn_out+"PROCESSED/FIGS/obs_phase_on_map.png")



if(0):
    # definitions for M2x and M2y are opposite to NEMO, but work:
    fig, [ax0, ax1] = plt.subplots(ncols=2)
    fes.dataset.M2x.plot(ax=ax0)
    nemo1.dataset.M2x.plot(ax=ax1)
    plt.show()


# Check that the phases are all aligned
list_of_datasets, list_of_titles = [obs, fes, nemo1, nemo2], ['obs', 'FES', 'GS1P1', 'GS1P7']
scatter_phase_on_map(list_of_datasets, list_of_titles)

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
ind5 = np.isnan(nemo1.dataset.M2x.values)
ind6 = np.isnan(nemo1.dataset.M2y.values)
ind7 = np.isnan(nemo2.dataset.M2x.values)
ind8 = np.isnan(nemo2.dataset.M2y.values)
#ind9 = np.isnan(tdiss.dataset.M2x.values)
#ind10 = np.isnan(tdiss.dataset.M2y.values)
I = ~ind1 * ~ind2 * ~ind3 * ~ind4 * ~ind5 * ~ind6 * ~ind7 * ~ind8 #* ~ind8 * ~ind10

#obs.dataset.M2x[~I] = np.nan
#obs.dataset.M2y[~I] = np.nan
#fes.dataset.M2x[~I] = np.nan
#fes.dataset.M2y[~I] = np.nan
#nemo1.dataset.M2x[~I] = np.nan
#nemo1.dataset.M2y[~I] = np.nan
#gs1p2.dataset.M2x[~I] = np.nan
#gs1p2.dataset.M2y[~I] = np.nan
obs.dataset = obs.dataset.isel(id_dim=I)
fes.dataset = fes.dataset.isel(id_dim=I)
nemo1.dataset = nemo1.dataset.isel(id_dim=I)
nemo2.dataset = nemo2.dataset.isel(id_dim=I)
#tdiss.dataset = tdiss.dataset.where(I)

print(ind1.flatten().sum())
print(ind2.flatten().sum())
print(ind3.flatten().sum())
print(ind4.flatten().sum())
print(ind5.flatten().sum())
print(ind6.flatten().sum())
print(ind7.flatten().sum())
print(ind8.flatten().sum())
#print(ind9.flatten().sum())
#print(ind10.flatten().sum())
print((~I).flatten().sum())


#%%  Plot distributions of depth at observation locations
plt.close('all')
plt.figure()
plt.plot( np.sort(np.log10(nemo1.dataset.bathymetry.values)) )
#plt.plot([0,500],[np.log10(2000), np.log10(2000)], 'm', label="2000m")
plt.plot([0,500],[np.log10(200), np.log10(200)], 'r', label="200m")
#plt.plot([0,500],[np.log10(50), np.log10(50)], 'g', label="50m")
#plt.plot([0,500],[np.log10(15), np.log10(15)], 'k', label="15m")
plt.xlabel('count')
plt.ylabel('log10(depth)')
plt.title("distribution of depths at observation sites")
plt.legend()
plt.savefig(config.dn_out+"PROCESSED/FIGS/dist_bathy.png")

obs, nemo1 = align_datasets(obs, nemo1)
obs, nemo2 = align_datasets(obs, nemo2)
obs, fes = align_datasets(obs, fes)
#obs, tdiss = align_datasets(obs, tdiss)

# define colors
c_ZPS_TIDE = 'lightblue' #[46/256., 153/256., 195/256.]
c_MES_TIDE = 'green'
c_FES = 'y'
c_obs = 'blue'


#%% Plot eaat coast of USA
#plot_east_coast_usa()

#%% ## Plot complex harmonic errors
#plot_all_complex_harmonic_errors()

#%% Plot the set of vector differences between model and obs for each constituent, coloured by depth
#plot_all_complex_harmonic_stick_errors()

#%% ## Plot amplitude and phase errors (deg)
#plot_all_amp_pha_errors()

#%% Compute Taylor diagram stats. One per constituent + depth class
#plot_all_taylor_tides()

#%% Compute Taylor diagrams. Overlay on deep and shallow plots as error trees
plot_overlay_taylor_tides()

#%% Attempt to do Taylor Tide with cloud of all data points
#plot_cloud()

#%% Make a table of output
#make_table()