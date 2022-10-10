import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from plotParam import *


## This code is to show the Distribution of colours in the mock catalogures.

fpath = './Data_Output/'
fname = 'galcatcat_pinocchio_mysel_test9.hdf5'#'galaxy_catalogue_pinocchio.hdf5'

fdata = np.loadtxt(fpath + fname + '.txt')
index_z = 2
## it can not be 0.

z_name = [0.0, 0.1, 0.2, 0.4, 0.6]
mag_name = [[17.75, 17.25], [18.75, 18.25], [19.75, 19.25], \
            [20.75, 20.25], [21.75, 21.25], [22.75, 22.25]]


#-----------Selection------------#
col = fdata[:, 3]
z_cos = fdata[:, 9]


z_cos_min = np.where(z_cos < z_name[index_z - 1])
z_cos_max = np.where(z_cos > z_name[index_z])

abs_mag = fdata[:, 0]
abs_mag_1775 = np.where(abs_mag < -17.75)
abs_mag_1725 = np.where(abs_mag > -17.25)
abs_mag_1875 = np.where(abs_mag < -18.75)
abs_mag_1825 = np.where(abs_mag > -18.25)
abs_mag_1975 = np.where(abs_mag < -19.75)
abs_mag_1925 = np.where(abs_mag > -19.25)
abs_mag_2075 = np.where(abs_mag < -20.75)
abs_mag_2025 = np.where(abs_mag > -20.25)
abs_mag_2175 = np.where(abs_mag < -21.75)
abs_mag_2125 = np.where(abs_mag > -21.25)
abs_mag_2275 = np.where(abs_mag < -22.75)
abs_mag_2225 = np.where(abs_mag > -22.25)


#------------Counts--------------#
nbin = 100
dcol = np.linspace(min(col), max(col), endpoint = True, num = nbin+1)

col_1 = np.delete(col, np.hstack((z_cos_min, z_cos_max, abs_mag_1775, abs_mag_1725)))
col_2 = np.delete(col, np.hstack((z_cos_min, z_cos_max, abs_mag_1875, abs_mag_1825)))
col_3 = np.delete(col, np.hstack((z_cos_min, z_cos_max, abs_mag_1975, abs_mag_1925)))
col_4 = np.delete(col, np.hstack((z_cos_min, z_cos_max, abs_mag_2075, abs_mag_2025)))
col_5 = np.delete(col, np.hstack((z_cos_min, z_cos_max, abs_mag_2175, abs_mag_2125)))
col_6 = np.delete(col, np.hstack((z_cos_min, z_cos_max, abs_mag_2275, abs_mag_2225)))


N1 = np.zeros(nbin); N2 = np.zeros(nbin); N3 = np.zeros(nbin)
N4 = np.zeros(nbin); N5 = np.zeros(nbin); N6 = np.zeros(nbin)
for i in range(nbin):
    N1[i] = np.count_nonzero(col_1 <= dcol[i+1])
    N2[i] = np.count_nonzero(col_2 <= dcol[i+1])
    N3[i] = np.count_nonzero(col_3 <= dcol[i+1])
    N4[i] = np.count_nonzero(col_4 <= dcol[i+1])
    N5[i] = np.count_nonzero(col_5 <= dcol[i+1])
    N6[i] = np.count_nonzero(col_6 <= dcol[i+1])

dN1 = np.zeros(nbin); dN2 = np.zeros(nbin); dN3 = np.zeros(nbin)
dN4 = np.zeros(nbin); dN5 = np.zeros(nbin); dN6 = np.zeros(nbin)
for i in range(nbin-1):
    dN1[i] = N1[i+1] - N1[i]
    dN2[i] = N2[i+1] - N2[i]
    dN3[i] = N3[i+1] - N3[i]
    dN4[i] = N4[i+1] - N4[i]
    dN5[i] = N5[i+1] - N5[i]
    dN6[i] = N6[i+1] - N6[i]


#-------------Plot---------------#
fig = plt.figure(1, (36, 20))
plt.plot(dcol[1:], dN1/180, label = r'$-17.25<^{0.1}M_r<-17.25$')
plt.plot(dcol[1:], dN2/180, label = r'$-18.25<^{0.1}M_r<-18.25$')
plt.plot(dcol[1:], dN3/180, label = r'$-19.25<^{0.1}M_r<-19.25$')
plt.plot(dcol[1:], dN4/180, label = r'$-20.25<^{0.1}M_r<-20.25$')
plt.plot(dcol[1:], dN5/180, label = r'$-21.25<^{0.1}M_r<-21.25$')
plt.plot(dcol[1:], dN6/180, label = r'$-22.25<^{0.1}M_r<-22.25$')


plt.xlabel(r'$^{0.1}(g-r)$', fontsize = 15)
plt.xlim([-0.2, 1.4])
plt.xticks([-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4], ('-0.2', '0.0', '0.2', '0.4', '0.6', '0.8', '1.0', '1.2'), fontsize = 15)
plt.ylabel(r'$N\,(deg^{-2})$', fontsize = 15)
#plt.ylim([0,5000])
#plt.yticks([0, 1000, 2000, 3000, 4000, 5000], ('0', '1000', '2000', '3000', '4000', '5000'), fontsize = 15)
#plt.ylim([0, 60000])
#plt.yticks([0, 10000, 20000, 30000, 40000, 50000, 60000], ('0', r'$1\times10^{4}$', r'$2\times10^{4}$', r'$3\times10^{4}$', r'$4\times10^{4}$', r'$5\times10^{4}$', r'$6\times10^{4}$'), fontsize = 13)

plt.text(1.0, 26000, '0.0 < z < '+str(z_name[index_z]), ha = 'right', ma = 'right', va = 'top', fontsize = 12) 
plt.legend(loc = 'upper left', fontsize = 10, frameon = False)
plt.tight_layout()
plt.savefig(fpath + fname + '_z_' + str(z_name[index_z]) + '_Ncol.png')
plt.show()


