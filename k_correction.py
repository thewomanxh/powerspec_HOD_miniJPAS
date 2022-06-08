#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import Func_Cosmo as Cfunc
from plotParam import *
#Copyright (c) Xiu-hui Tan
## This code is to build k correction for JPAS.

fpath = './Data_Input/JPAS/'
fcsv = 'minijpas_pdr201912_cat_wAbsMag_v1.csv'
fout = './Data_Output/'

h = 0.6774
def calc_kcor(izz, rapp, iabs):
    Dl = Cfunc.Distfunc(izz, unit = 'Mpc', ih = h)
    k_04 = rapp - 5*np.log10(Dl) - 25. - iabs + 5*np.log10(h)
    return k_04

#a=1.4+5.*np.log10(2.64/10)
#a=-1.46+5.*np.log10(10/2.64)

fdata = pd.read_csv(fpath + fcsv)
fdata_array = np.array(fdata)
iband_abs = fdata.get('i_abs') #- 5*np.log10(h)
rband_app = fdata.get('mag_r_auto')
g_r = fdata.get('g_abs') - fdata.get('r_abs')
photoz = fdata.get('PHOTOZ')

## delete the stars with CLASS_STAR <= 0.5
stars = fdata.get('CLASS_STAR')
stars_index = (stars[stars.values<=0.5]).index
iband_abs = iband_abs.drop(stars_index)
rband_app = rband_app.drop(stars_index)
g_r = g_r.drop(stars_index)
photoz = photoz.drop(stars_index)

## delete z=0
z0_index = (photoz[photoz.values==0]).index
#z02 = (photoz[(photoz.values<0.2)]).index
#z06 = (photoz[(photoz.values>0.6)]).index
#z0_index = z02.append(z06)
iband_abs = iband_abs.drop(z0_index)
rband_app = rband_app.drop(z0_index)
g_r = g_r.drop(z0_index)
zz = photoz.drop(z0_index)

## bin g_r range
#g_r_range = [-100, -0.18, 0.35, 0.52, 0.69, 0.86, 1.03, 100.]
g_r_range = np.linspace(min(g_r), max(g_r), num = 8)
g_r_bin_index = [[] for i in range(7)]

plt.figure(1, (18, 16))
for i in range(7):
    aa = (g_r[(g_r > g_r_range[i]) & (g_r <= g_r_range[i+1])]).index
    g_r_bin_index[i].append(aa)
    index = np.array(g_r_bin_index[i][0])
    k0 = calc_kcor(np.array(zz[index]), np.array(rband_app[index]), np.array(iband_abs[index]))
    #kcor = np.array(np.column_stack([np.array(photoz[index]), k0]))
    #print(zz, k0)

    ## plot points.
    plt.plot(zz[index], k0, '.')
    plt.xlim([0, 1.5])
    plt.xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4], (r'$0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$', r'$1.0$', r'$1.2$', r'$1.4$'), fontsize = 15)
    plt.xlabel(r'$z$', fontsize = 17)
    plt.ylabel(r'$k(z)$', fontsize = 17)
    plt.tight_layout()
    plt.show()
   
    ## assign polynomial
    z1 = np.polyfit(zz[index], k0, 4)
    p1 = np.poly1d(z1)
    print(i, '\n', p1)
    
    #y = lambda x: z1[0]*x**4 + z1[1]*x**3 + z1[2]*x**2 + z1[3]*x + z1[4]
    xx = np.arange(0, 1.5, 0.02)
    plt.plot(xx, p1(xx), '-')
    plt.xlim([0, 1.5])
    plt.xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4], (r'$0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$', r'$1.0$', r'$1.2$', r'$1.4$'), fontsize = 15)
    plt.xlabel(r'$z$', fontsize = 17)
    plt.ylabel(r'$k(z)$', fontsize = 17)
    plt.tight_layout()
    plt.show()



