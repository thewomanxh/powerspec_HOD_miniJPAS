import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import interp1d,interp2d
from plotParam import *


## This code is to show the dNdz in the mock catalogures.

fpath = './Data_Output/'
fname = 'galcatcat_pinocchio_mysel_test9.hdf5'#'galaxy_catalogue_pinocchio.hdf5'

fdata = np.loadtxt(fpath + fname + '.txt')
bin_lum = ['-25.2', '-23.0', '-22.0', '-21.0', '-20.0', '-17.5']
mag_name = [-25.2, -23.0, -22.0, -21.0, -20.0, -17.5]
ibin_index = 1
zbins = np.linspace(0.0001, 2., num = 100)
zbins_mid = np.zeros(len(zbins) - 1)
dndz = np.zeros(len(zbins) - 1)
dndz_deltazz = np.zeros(len(zbins) - 1)

abs_mag = fdata[:, 0]
abs_mag_del = np.where(abs_mag > mag_name[ibin_index + 1])

z_cos = fdata[:, 9]
#print(z_cos);exit()
z_dist = np.delete(z_cos, abs_mag_del)

#for i in range(len(zbins) - 1):
#    zbins_mid[i] = (zbins[i] + zbins[i+1])/2.
#    deltazz = (zbins[i+1] - zbins[i])
#    dndz[i] = (((zbins[i]<z_dist) & (z_dist<zbins[i+1])).sum())
#    dndz_deltazz[i] = dndz[i]/deltazz
#winres = dndz_deltazz/sum(dndz)


#-------------Plot---------------#
fig = plt.figure(1, (36, 20))
plt.hist(z_dist, bins = 100, density = True)
plt.ylabel(r'$\mathrm{d}N\mathrm{d}z$', fontsize = 15)
#plt.ylabel(r'$\mathrm{W\,function}$', fontsize = 15)
plt.xlabel(r'$z$', fontsize = 15)
#plt.xlim([0.2, 0.6])
#plt.xticks([0.2, 0.3, 0.4, 0.5, 0.6], ('0.2', '0.3', '0.4', '0.5', '0.6'), fontsize = 12)
plt.legend(loc = 'upper left', fontsize = 10, frameon = False)
plt.tight_layout()
plt.savefig('./Data_Output/miniJPAS_axiliary/mock_' + bin_lum[0] + '_' + bin_lum[ibin_index+1] + '_dndz.png')
plt.show()


