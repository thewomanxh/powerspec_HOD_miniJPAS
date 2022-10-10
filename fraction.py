#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Copyright (c) Xiu-hui Tan
## This code is to build the fraction of blue galaxies and central galaxies for AMEGO, which is in name of miniJPAS.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fpath = '../Data_Input/AMEGO/'
f_cenblue = 'galaxies_central_blue.csv'
f_cenred = 'galaxies_central_red.csv'
f_satblue = 'galaxies_satelite_blue.csv'
f_satred = 'galaxies_satelite_red.csv'
fout = './Data_Output/'

def read_data(fname, Sat = False, Out = False):
    fdata = np.array(pd.read_csv(fname))
    ra = fdata[:, 0]
    dec = fdata[:, 1]
    mag_auto_r = fdata[:, 2]
    photonz = fdata[:, 3]
    if Sat:
        id_assoc = None
        prob_assoc = fdata[:, 4]
    else:
        id_assoc = fdata[:, 4]
        prob_assoc = fdata[:, 5]
    if Out:
        return fdata[:, Out]
    else:
        return ra, dec, mag_auto_r, photonz, id_assoc, prob_assoc

ra, dec, mag_satblue, z_satblue, x, prob_satblue = read_data(fpath + f_satblue, Sat = True)
ra, dec, mag_cenblue, z_cenblue, id_cenblue, prob_cenblue = read_data(fpath + f_cenblue)
ra, dec, mag_cenred, z_cenred, id_cenred, prob_cenred = read_data(fpath + f_cenred)
ra, dec, mag_satred, z_satred, x, prob_satred = read_data(fpath + f_satred, Sat = True)
print(mag_satblue, z_satblue, len(prob_satblue) )
print(mag_cenblue, z_cenblue, prob_cenblue)

ax = plt.subplot(111, projection = '3d')
ax.scatter(mag_satblue, z_satblue, prob_satblue)
ax.set_zlabel('prob')
ax.set_xlabel('mag')
ax.set_ylabel('redshift')
plt.show()





