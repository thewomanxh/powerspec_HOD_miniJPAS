import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from plotParam import *
from compute_HOD import dndm_func


## This code is to show the dNdz in the mock catalogures.

fpath = './Data_Output/'
fname = 'galcatcat_pinocchio_mysel_test9.hdf5'#'galaxy_catalogue_pinocchio.hdf5'

fdata = np.loadtxt(fpath + fname + '.txt')
bin_lum = ['-25.2', '-23.0', '-22.0', '-21.0', '-20.0', '-17.5']
mag_name = [-25.2, -23.0, -22.0, -21.0, -20.0, -17.5]
zz_fix = 0.4
zcos = fdata[:, 9]
zz05_del = np.where(zcos>0.5)
zz03_del = np.where(zcos<0.3)
halo_mass = np.delete(fdata[:, 6], np.hstack((zz03_del, zz05_del)))
abs_mag = np.delete(fdata[:, 0], np.hstack((zz03_del, zz05_del)))
dndm = dndm_func(halo_mass, zz_fix)*halo_mass

#-------------Plot---------------#
fig = plt.figure(1, (36, 20))
plt.plot(abs_mag, dndm)#, lable = 'z ='+str(zz_fix))
plt.ylabel(r'$\phi\,h^3\mathrm{Mpc^{-3}mag^{-1}}$', fontsize = 15)
plt.xlabel(r'$^{0.4}M_i$', fontsize = 15)
#plt.xlim([-25., -17.])
#plt.legend(loc = 'upper left', fontsize = 10, frameon = False)
#plt.tight_layout()
#plt.savefig('./Data_Output/miniJPAS_axiliary/mock_z_' + str(zz_fix) + '_luminosity.png')
plt.show()


