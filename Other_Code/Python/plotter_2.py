#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 22:02:25 2017

@author: ellismle
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle



plot_index = 'end'



colours = ['r', 'g', 'b']
Plot_XYZ = 'xyz'

filepath = '../../Data/'

def p(y):
    plt.close()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x,y)
    ax_filler(ax)
    ax.set_xlabel('X [nm]')
    fig.show()


def read_data(path = filepath, files = ['Ux', 'Uy', 'Uz', 'Jmx', 'Jmy', 'Jmz', 'params_vec', 'params_scalar']):   
    if path[-1] != '/':
        path += '/'
    filepaths = [path + i for i in files]
    data = [pd.read_csv(i) for i in filepaths]    
    return data
ux_orig, uy_orig, uz_orig, Jmx_orig, Jmy_orig, Jmz_orig, params_v, params_s = read_data()
ux = ux_orig.drop('dt', 1)
uy = uy_orig.drop('dt', 1)
uz = uz_orig.drop('dt', 1)
Jmx = Jmx_orig.drop('dt', 1)
Jmy = Jmy_orig.drop('dt', 1)
Jmz = Jmz_orig.drop('dt', 1)

def plot_index_convert(ind = plot_index, U = uz):
    if 'end' in ind:
        ind = len(U)-1
    return ind

def ax_filler(ax, intensity=0.15):
    color_intensities = 1 - (beta/np.max(beta)* intensity)
    for i in range(len(Uz)):
        if beta[i] > 0.01:
            ax.add_artist(Rectangle((i*dx*1e9,-1e20),dx*1e9,2e20, facecolor=str(color_intensities[i]), edgecolor="none"))
    return ax


time_steps = ['dt']
x = params_v['x']*1e9
beta = params_v['Beta']
beta_prime = params_v['Beta_Prime']
D = params_v['D']
lambda_sf = params_v['Lambda_Sf']
lambda_J = params_v['Lambda_J']
lambda_phi = params_v['Lambda_Phi']
Mx = params_v['Mx']
My = params_v['My']
Mz = params_v['Mz']
dx = params_s['dx'][0]
dt = params_s['dt'][0]
L = params_s['L'][0]
T = params_s['T'][0]
DL = params_s['Dl'][0]
j_e = params_s['je'][0]
Mx_peaks = params_s['Mx']
My_peaks = params_s['My']
Mz_peaks = params_s['Mz']
lambda_sf_peaks = params_s['Lambda_Sf']
lambda_J_peaks = params_s['Lambda_J']
lambda_phi_peaks = params_s['Lambda_Phi']
beta_peaks = params_s['Beta']
beta_prime_peaks = params_s['Beta_Prime']
D_peaks = params_s['D']

plot_index = plot_index_convert()

Uz = uz.loc[plot_index]
Uy = uy.loc[plot_index]
Ux = ux.loc[plot_index]

jmx = Jmx.loc[plot_index]
jmy = Jmy.loc[plot_index]
jmz = Jmz.loc[plot_index]

def peak_finders(u):
    mult = 3
    try:
        max_x = x[int(np.argmax(u[u > 1]))]
        max_x  = [max_x + (mult*lambda_sf_peaks[0])*1e9, max_x - (mult*lambda_sf_peaks[0])*1e9]
        print(max_x)
    except ValueError:
        max_x = [False]
    try: 
        min_x =x[int( np.argmin(u[u < -1]))]
        min_x  = [min_x +(mult*lambda_sf_peaks[0])*1e9, min_x - (mult*lambda_sf_peaks[0])*1e9]
    except ValueError:
        min_x = [False]
    tot_peak_locs = max_x + min_x
    return [np.min(tot_peak_locs), np.max(tot_peak_locs)]


def plot_setup():
    fig = plt.figure(facecolor='white', figsize=(16, 8))
    ax = plt.subplot2grid((6,2), (0,0), rowspan=5)
    ax_jm = plt.subplot2grid((6,2), (0,1), rowspan=5)
    ax_params = plt.subplot2grid((6,2), (5,0), colspan=2)
    params = [beta_peaks, beta_prime_peaks, D_peaks, lambda_sf_peaks, lambda_phi_peaks, lambda_J_peaks, Mx_peaks, My_peaks, Mz_peaks]
    ax_params.axis('off')
    param_data = [['Non-Magnet']+["%.2g"%i[0] for i in params], ['F1']+["%.2g"%i[1] for i in params], ['F2']+["%.2g"%i[2] for i in params] ]
    headers = ['',r'$\beta$', r"$\beta$'",r"$D$", r"$\lambda_{sf}$", r"$\lambda_{\phi}$", r"$\lambda_{J}$", r"$M_x$", r"$M_y$", r"$M_z$"]
    table = ax_params.table(cellText = param_data, loc='center', colLabels=headers, cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(15)
    table.scale(1, 2.5)
    
    plt.suptitle(r"DL = %.2gnm    j$_e$ = %.2g Am$^{-2}$     T = %.2gfs"%(DL*1e9, j_e, T*1e15), fontsize = 16)    
    
    ax_filler(ax)
    ax.set_title("Spin Accumuation -Z", fontsize=14)
    ax.set_xlabel(r"X [nm]")
    ax_jm.set_xlabel(r"X [nm]")
    #ax.set_ylabel(r"Uz [$C m ^{-3}$]")
    #ax_jm.set_ylabel(r"$j_{m_{z}}$ [$A m^{-2}$]")
    ax.set_xlim([0,L*1e9])
    
    ax_filler(ax_jm)
    ax_jm.set_title("Spin Current -Z", fontsize=14)
    xlims = peak_finders(Uz)
    
    
    ax.set_xlim(xlims)
    ax_jm.set_xlim(xlims)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return fig, ax, ax_jm, ax_params

fig, Uax, Jmax, ax_params = plot_setup()

if 'x' in Plot_XYZ.lower():
    Uax.plot(x, Ux, colours[0])
    Jmax.plot(x, jmx, colours[0])
if 'y' in Plot_XYZ.lower():
    Uax.plot(x, Uy, colours[1])
    Jmax.plot(x, jmy, colours[1])
if 'z' in Plot_XYZ:
    Uax.plot(x, Uz, colours[2])
    Jmax.plot(x, jmz, colours[2])

def test_params():
    def get_variable_name(*variable):
        if len(variable) != 1:
            raise Exception('len of variables inputed must be 1')
        try:
            return [k for k, v in locals().items() if v is variable[0]][0]
        except:
            return [k for k, v in globals().items() if v is variable[0]][0]

    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    ax.set_xlabel("X [nm]")
    test_params = [beta, beta_prime, D, lambda_sf, lambda_phi, lambda_J, Mx, My, Mz]
    for i in test_params:
        ax.clear()
        ax.plot(x, i)
        ax.set_ylabel(get_variable_name(i))
        plt.pause(3)

plt.show()

#fig.savefig('/home/ellismle/Documents/Mphys_Project/Graphs/Diffuse_interfaces/Spin_Acc_Z_3nm.png', dpi=900, format='png')
#fig1.savefig('/home/ellismle/Documents/Mphys_Project/Graphs/Diffuse_interfaces/Spin_Current_Z_3nm.png', dpi=900, format='png')
