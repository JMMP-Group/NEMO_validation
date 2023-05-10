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


def amp_pha_from_re_im(creal, cimag):
    """
    Example usage: amp,pha = amp_pha_from_re_im(ds.M2x, ds.M2y)
    """
    if type(creal) == xr.DataArray:
        amp = xr.zeros_like(creal)
        pha = xr.zeros_like(creal)
    cc=creal+cimag*1j
    amp=np.abs(cc)
    pha=np.angle(cc)*180/np.pi
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
    for EXP in ["GS1P1", "GS1P2", "FES2014"]:
        for constit in constit_list:
            if EXP == "GS1P1": tg_mod = gs1p1
            if EXP == "GS1P2": tg_mod = gs1p2
            if EXP == "FES2014": tg_mod = fes

            # separate observations by depth
            try:
                ind_deep = tg_mod.dataset.bathymetry.values > 200
            except:
                ind_deep = gs1p1.dataset.bathymetry.values > 200
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
    for EXP in ["GS1P1", "GS1P2", "FES2014"]:
        for constit in constit_list:
            if EXP == "GS1P1": tg_mod = gs1p1
            if EXP == "GS1P2": tg_mod = gs1p2
            if EXP == "FES2014": tg_mod = fes

            # separate observations by depth
            try:
                ind_deep = tg_mod.dataset.bathymetry.values > 200
            except:
                ind_deep = gs1p1.dataset.bathymetry.values > 200
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
    for EXP in ["GS1P1", "GS1P2", "FES2014"]:
        for constit in constit_list:
            if EXP == "GS1P1": tg_mod = gs1p1
            if EXP == "GS1P2": tg_mod = gs1p2
            if EXP == "FES2014": tg_mod = fes
            print(f"Plot {constit} amplitude and phase errors (deg): {EXP}")

            # separate observations by depth
            try:
                ind_deep = tg_mod.dataset.bathymetry.values > 200
            except:
                ind_deep = gs1p1.dataset.bathymetry.values > 200
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
                II = gs1p1.dataset.bathymetry.values > 200
            elif subset == 'shal':
                II = gs1p1.dataset.bathymetry.values <= 200
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
            if count == 0:
                R = np.array([1])  # R for obs
                rms_amp = np.array([rms_abs_amp(z1obs, z2obs)])
                rms_err = np.array([0])  # err for obs
            else:
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

            # GS1P1
            del z1mod, z2mod
            z1mod, z2mod = gs1p1.dataset[constit + 'x'][II], gs1p1.dataset[constit + 'y'][II]
            # z1obs_new, z1mod_new = tganalysis.match_missing_values(obs.dataset.M2x, tg_mx2.dataset.M2x)
            # z2obs_new, z2mod_new = tganalysis.match_missing_values(obs.dataset.M2y, tg_mx2.dataset.M2y)

            R = np.hstack((R, pearson_correl_coef(z1obs, z2obs, z1mod, z2mod)))
            rms_amp = np.hstack((rms_amp, rms_abs_amp(z1mod, z2mod)))
            rms_err = np.hstack((rms_err, rms_abs_error(z1obs, z2obs, z1mod, z2mod)))

            # GS1P2
            del z1mod, z2mod
            try:
                z1mod, z2mod = gs1p2.dataset[constit + 'x'][II], gs1p2.dataset[constit + 'y'][II]
            except:
                z1mod = np.nan
                z2mod = np.nan
            # z1obs_new, z1mod_new = tganalysis.match_missing_values(obs.dataset.M2x, tg_mv2.dataset.M2x)
            # z2obs_new, z2mod_new = tganalysis.match_missing_values(obs.dataset.M2y, tg_mv2.dataset.M2y)
            R = np.hstack((R, pearson_correl_coef(z1obs, z2obs, z1mod, z2mod)))
            rms_amp = np.hstack((rms_amp, rms_abs_amp(z1mod, z2mod)))
            rms_err = np.hstack((rms_err, rms_abs_error(z1obs, z2obs, z1mod, z2mod)))

            count = count + 1

        label = ['obs:s', 'fes:s', 'gs1p1:s', 'gs1p2:s',
                 'obs:d', 'fes:d', 'gs1p1:d', 'gs1p2:d']

        print(f"R= {[format(R[i], '.2f') for i in range(len(R))]}")
        print(f"rms_amp= {[format(rms_amp[i], '.2f') for i in range(len(R))]}")
        print(f"rms_err= {[format(rms_err[i], '.2f') for i in range(len(R))]}")
        print(label)

        ## Check cosine rule consistency
        A = rms_amp
        C = rms_err
        costheta = R

        B = rms_amp[0]
        for i in range(0, 4):
            print(
                f"{label[i]}: sqrt(A^2+B^2-2ABcos(theta))={np.sqrt(A[i] ** 2 + B ** 2 - 2 * A[i] * B * costheta[i])}. C={C[i]}")
        B = rms_amp[5]
        for i in range(4, 8):
            print(
                f"{label[i]}: sqrt(A^2+B^2-2ABcos(theta))={np.sqrt(A[i] ** 2 + B ** 2 - 2 * A[i] * B * costheta[i])}. C={C[i]}")

        # Create TaylorTide plot template
        tt = TaylorTide(
            r_obs=rms_amp[0],
            rms_amp_max=0.7,
            rms_amp_contours=[0.2, 0.4, 0.6],
            rms_err_contours=[0.2, 0.4, 0.6],
            cos_theta_lines=[0.3, 0.6, 0.9],
        )
        # Add data to axes
        tt.ax.scatter(rms_amp[1:4] * R[1:4], rms_amp[1:4] * np.sqrt(1 - R[1:4] ** 2), s=10, c=['r', 'k', 'g'])
        # manual legend
        colors = ['red', 'black', 'green']
        lines = [Line2D([0], [0], color=c, linewidth=3, linestyle='dotted') for c in colors]
        labels = ["FES2014", "GS1p1", "GS1p2"]
        plt.legend(lines, labels)

        plt.title(constit + ':shallow')
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
        tt.ax.scatter(rms_amp[5::] * R[5::], rms_amp[5::] * np.sqrt(1 - R[5::] ** 2), s=10, c=['r', 'k', 'g'])
        # manual legend
        colors = ['red', 'black', 'green']
        lines = [Line2D([0], [0], color=c, linewidth=3, linestyle='dotted') for c in colors]
        labels = ["FES2014", "GS1p1", "GS1p2"]
        plt.legend(lines, labels)

        plt.title(constit + ':deep')
        plt.savefig(config.dn_out + "PROCESSED/FIGS/Taylor_" + constit + "_deep.png")

def plot_overlay_taylor_tides():
    # Loop over harmonic species.
    # Plot Taylor Tide diag of model and obs for each harmonic. Overlay on two (deep/shallow) plots as trees
    for subset in ['shal', 'deep']:
      for constit_family_list in [["M2", "S2", "N2", "K2"], ["O1", "Q1", "P1"]]:
        if "M2" in constit_family_list: family_str = "semi-diurnal"
        if "O1" in constit_family_list: family_str = "diurnal"
        R = np.zeros((4, len(constit_family_list)))
        rms_amp = np.zeros((4, len(constit_family_list)))
        rms_err = np.zeros((4, len(constit_family_list)))
        label = {}
        for count, constit in enumerate(constit_family_list):
            try:
                del II, z1obs, z2obs, z1mod, z2mod
            except:
                pass
            if subset == 'deep':
                # separate observations by depth
                II = gs1p1.dataset.bathymetry.values > 200
            elif subset == 'shal':
                II = gs1p1.dataset.bathymetry.values <= 200
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


            # GS1P1
            try:
                del z1mod, z2mod
            except:
                pass
            z1mod, z2mod = gs1p1.dataset[constit + 'x'][II], gs1p1.dataset[constit + 'y'][II]

            R[1,count] = pearson_correl_coef(z1obs, z2obs, z1mod, z2mod)
            rms_amp[1,count] = rms_abs_amp(z1mod, z2mod)
            rms_err[1,count] = rms_abs_error(z1obs, z2obs, z1mod, z2mod)
            label[1,count] = 'gs1p1:'+constit



            # GS1P2
            del z1mod, z2mod
            #try:
            z1mod, z2mod = gs1p2.dataset[constit + 'x'][II], gs1p2.dataset[constit + 'y'][II]
            #except:
            #    z1mod = np.nan
            #    z2mod = np.nan

            R[2,count] = pearson_correl_coef(z1obs, z2obs, z1mod, z2mod)
            rms_amp[2,count] = rms_abs_amp(z1mod, z2mod)
            rms_err[2,count] = rms_abs_error(z1obs, z2obs, z1mod, z2mod)
            label[2,count] = 'gs1p2:'+constit


            # FES
            del z1mod, z2mod
            z1mod, z2mod = fes.dataset[constit + 'x'][II], fes.dataset[constit + 'y'][II]

            R[3,count] = pearson_correl_coef(z1obs, z2obs, z1mod, z2mod)
            rms_amp[3,count] = rms_abs_amp(z1mod, z2mod)
            rms_err[3,count] = rms_abs_error(z1obs, z2obs, z1mod, z2mod)
            label[3,count] = 'fes:'+constit

        print(subset)
        #print(f"R= {[format(R[i], '.2f') for i in range(len(R))]}")
        #print(f"rms_amp= {[format(rms_amp[i], '.2f') for i in range(len(R))]}")
        #print(f"rms_err= {[format(rms_err[i], '.2f') for i in range(len(R))]}")
        print(label)

        ## Check cosine rule consistency
        for i in range(len(constit_family_list)):
            B = rms_amp[0,i]
            A = rms_amp[1:4,i]
            C = rms_err[1:4,i]
            costheta = R[1:4,i]

            for j in range(3): # model runs: fes, gs1p1, gs1p2
                print(
                    f"{label[1+j,i]}: sqrt(A^2+B^2-2ABcos(theta))={np.sqrt(A[j] ** 2 + B ** 2 - 2 * A[j] * B * costheta[j])}. C={C[j]}")
            del B, A, C, costheta




        # Create TaylorTide plot template
        tt = TaylorTide(
            r_obs=rms_amp[0,0],
            rms_amp_max=0.7,
            rms_amp_contours=[0.2, 0.4, 0.6],
            rms_err_contours=[0.2, 0.4, 0.6],
            cos_theta_lines=[0.3, 0.6, 0.9],
            err_contour_flag=False,
        )

        ## Loop over constituents to create plot
        for i in range(len(constit_family_list)):
            # Add data to axes
            tt.ax.scatter(rms_amp[0:4,i] * R[0:4,i],
                          rms_amp[0:4,i] * np.sqrt(1 - R[0:4,i] ** 2),
                          s=20, c=['b', 'k', 'g', 'r'])
            # Add vectors between points and obs
            tt.ax.plot([np.repeat(rms_amp[0,i],3), rms_amp[1:4,i] * R[1:4,i]],
                       [np.zeros(3), rms_amp[1:4,i] * np.sqrt(1 - R[1:4,i] ** 2)])
            tt.ax.lines[-3].set_color('k')
            tt.ax.lines[-2].set_color('g')
            tt.ax.lines[-1].set_color('r')

            tt.ax.text( rms_amp[0,i], -0.025, constit_family_list[i], rotation=0, color='b')


        # manual legend
        colors = ['black', 'green', 'red', 'blue']
        lines = [Line2D([], [], color=c, markersize=5, marker='o', linestyle='None') for c in colors]
        labels = ["GS1p1", "GS1p2", "FES2014", "obs"]
        plt.legend(lines, labels, loc='upper left')

        plt.title(subset + ":" + family_str)
        plt.savefig(config.dn_out + "PROCESSED/FIGS/Taylor_" + "_" + subset + "_" + family_str + "_tree.png")

def plot_cloud():
    # Taylor Tide with cloud of all data points
    #for constit in constit_list:
    constit = "M2"
    EXP = "GS1P1" #, "GS1P2", "FES2014"]:
    if EXP == "GS1P1": tg_mod = gs1p1
    if EXP == "GS1P2": tg_mod = gs1p2
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
            II = gs1p1.dataset.bathymetry.values > 200
        elif subset == 'shal':
            II = gs1p1.dataset.bathymetry.values <= 200
        else:
            print(f"Not expecting that {subset}")

        z1obs, z2obs = obs.dataset[constit + 'x'][II], obs.dataset[constit + 'y'][II]

        # Obs
        if count == 0:
            R = np.array([1])  # R for obs
            rms_amp = np.array([rms_abs_amp(z1obs, z2obs)])
            rms_amp_obs = np.array([rms_abs_amp(z1obs, z2obs)])
            rms_err = np.array([0])  # err for obs


        z1mod, z2mod = tg_mod.dataset[constit + 'x'][II], tg_mod.dataset[constit + 'y'][II]

        # model - averaged
        if count == 1:
            R = np.hstack((R, pearson_correl_coef(z1obs, z2obs, z1mod, z2mod)))
            rms_amp = np.hstack((rms_amp, rms_abs_amp(z1mod, z2mod)))
            rms_amp_obs = np.hstack((rms_amp_obs, rms_abs_amp(z1obs, z2obs)))
            rms_err = np.hstack((rms_err, rms_abs_error(z1obs, z2obs, z1mod, z2mod)))

        n_length = len(z1obs)
        for ii in range(n_length):
            R = np.hstack((R, pearson_correl_coef(z1obs[ii], z2obs[ii], z1mod[ii], z2mod[ii])))
            rms_amp = np.hstack((rms_amp, rms_abs_amp(z1mod[ii], z2mod[ii])))
            rms_amp_obs = np.hstack((rms_amp_obs, rms_abs_amp(z1obs[ii], z2obs[ii])))
            rms_err = np.hstack((rms_err, rms_abs_error(z1obs[ii], z2obs[ii], z1mod[ii], z2mod[ii])))




        label = ['obs:s', 'fes:s', 'gs1p1:s', 'gs1p2:s',
                 'obs:d', 'fes:d', 'gs1p1:d', 'gs1p2:d']

        print(f"R= {[format(R[i], '.2f') for i in range(len(R))]}")
        print(f"rms_amp= {[format(rms_amp[i], '.2f') for i in range(len(R))]}")
        print(f"rms_err= {[format(rms_err[i], '.2f') for i in range(len(R))]}")
        print(label)

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
            rms_amp_max=1.75,
            rms_amp_contours=[0.2, 0.4, 0.6, 1, 1.5],
            rms_err_contours=[0.2, 0.4, 0.6],
            cos_theta_lines=[0.3, 0.6, 0.9],
        )
        # Add data to axes
        tt.ax.scatter(rms_amp[2:n_length] * R[2:n_length], rms_amp[2:n_length] * np.sqrt(1 - R[2:n_length] ** 2), s=10, c='k')
        tt.ax.scatter(rms_amp[0] * R[0], rms_amp[0] * np.sqrt(1 - R[0] ** 2), s=20, c='g')
        tt.ax.scatter(rms_amp[1] * R[1], rms_amp[1] * np.sqrt(1 - R[1] ** 2), s=20, c='r')
        tt.ax.scatter(rms_amp[2] * R[2], rms_amp[2] * np.sqrt(1 - R[2] ** 2), s=20, c='m')

        for i in [0,1,2]:
            print(f"{rms_amp[i] * R[i], rms_amp[i] * np.sqrt(1 - R[i]**2)}")

        plt.title(EXP + " " + constit + ':dots')
        plt.savefig(config.dn_out + "PROCESSED/FIGS/Taylor_" + constit + "_" + EXP + "_shallow_SPECIAL.png")

def plot_east_coast_usa():
    ylims = [35, 45]  # [-80,80] #[45,60]
    xlims = [-76, -60]  # [-180,180] #[-60, -25]

    plt.close('all')
    #plt.figure()
    fig, axs = plt.subplots(2, 2, figsize=(10, 5), subplot_kw={'projection': ccrs.PlateCarree()})
    axs[0,0].coastlines()
    axs[0,0].scatter(obs.dataset.longitude, obs.dataset.latitude, c=obs.dataset['A'], s=20,
                     vmin=0, vmax=0.45, cmap='Spectral')
    axs[0,0].set_title('obs')

    axs[0,1].coastlines()
    axs[0,1].scatter(fes.dataset.longitude, fes.dataset.latitude, c=fes.dataset['A'], s=20,
                     vmin=0, vmax=0.45, cmap='Spectral')
    axs[0,1].set_title('FES')


    axs[1,0].coastlines()
    axs[1,0].scatter(gs1p1.dataset.longitude, gs1p1.dataset.latitude, c=gs1p1.dataset['A'], s=20,
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

# Load data as tidegauge objects
## Harmonise definitions: negate M2y (and phase) in NEMO - done in preprocessing.
################################
obs = coast.Tidegauge(dataset=xr.open_dataset(config.dn_out+"PROCESSED/obs_extracted.nc"))
obs.dataset['A'], obs.dataset['G'] = amp_pha_from_re_im(obs.dataset.M2x, obs.dataset.M2y)
obs.dataset['longitude'] = ((obs.dataset.longitude%360) + 180)%360 - 180  # set longitude : [-180,180]
obs.dataset['G'] = (obs.dataset.G+180)%360-180  # set phase: -180,180

fes = coast.Tidegauge(dataset=xr.open_dataset(config.dn_out+"PROCESSED/FES2014_extracted.nc"))
fes.dataset['A'], fes.dataset['G'] = amp_pha_from_re_im(fes.dataset.M2x, fes.dataset.M2y)
fes.dataset['G'] = (fes.dataset.G+180)%360-180  # set phase: -180,180

gs1p1 = coast.Tidegauge(dataset=xr.open_dataset(config.dn_out+"PROCESSED/GS1p1_tide_extracted.nc"))
gs1p1.dataset['A'], gs1p1.dataset['G'] = amp_pha_from_re_im(gs1p1.dataset.M2x, gs1p1.dataset.M2y)

gs1p1.dataset['G'] = (gs1p1.dataset.G+180)%360-180  # set phase: -180,180
gs1p1.dataset['M2y'] = +gs1p1.dataset.M2y
gs1p1.dataset = gs1p1.dataset.drop_dims(["nvertex", "z_dim"])  # drop unwanted dimensions and associated variables

gs1p2 = coast.Tidegauge(dataset=xr.open_dataset(config.dn_out+"PROCESSED/GS1p2_full_extracted.nc"))
#gs1p2.dataset = gs1p2.dataset.drop_vars(["G", "A"]) # Unlabelled amp/pha parts should not be in the file
gs1p2.dataset['A'], gs1p2.dataset['G'] = amp_pha_from_re_im(gs1p2.dataset.M2x, gs1p2.dataset.M2y)

gs1p2.dataset['G'] = -(gs1p2.dataset.G+180)%360-180  # set phase: -180,180
gs1p2.dataset['M2y'] = +gs1p2.dataset.M2y
gs1p2.dataset = gs1p2.dataset.drop_dims(["nvertex", "z_dim"])  # drop unwanted dimensions and associated variables



# Compare maps
if(0):
 # Plot utility for COAsT.tidegauge objects
 obs.plot_on_map_multiple([obs], color_var_str="G")
 #plt.savefig(config.dn_out+"PROCESSED/FIGS/obs_phase_on_map.png")

ylims = [-80,80] #[45,60]
xlims = [-180,180] #[-60, -25]

if(0):
    # definitions for M2x and M2y are opposite to NEMO, but work:
    fig, [ax0, ax1] = plt.subplots(ncols=2)
    fes.dataset.M2x.plot(ax=ax0)
    gs1p1.dataset.M2x.plot(ax=ax1)
    plt.show()

## Basic map plot of phases
# Phase alignment

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

#%% Plot eaat coast of USA
plot_east_coast_usa()

#%% ## Plot complex harmonic errors
#plot_all_complex_harmonic_errors()

#%% Plot the set of vector differences between model and obs for each constituent, coloured by depth
#plot_all_complex_harmonic_stick_errors()

#%% ## Plot amplitude and phase errors (deg)
#plot_all_amp_pha_errors()

#%% Compute Taylor diagram stats. One per constituent + depth class
#plot_all_taylor_tides()

#%% Compute Taylor diagrams. Overlay on deep and shallow plots as error trees
#plot_overlay_taylor_tides()

#%% Attempt to do Taylor Tide with cloud of all data points
#plot_cloud()