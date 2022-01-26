import numpy as np
import  matplotlib.pyplot as plt

## This code is to divide miniJPAS field into Nsamples subregions.


#---------Global parameters-------#
global bin_lum
bin_lum = ['-25.2', '-17.5', '-20.0', '-21.0', '-22.0', '-23.0']
bin_num = 0
Nsamples = 20


## Divide Nsamples samples.
binfile = np.loadtxt('./Data_Output/miniJPAS_redshift/miniJPAS_zz_nostar_iband_selection_' + bin_lum[0] + '_' + bin_lum[bin_num+1] + '.txt')
coor = binfile[:, 1:3]
coor1 = binfile[:, 1]
coor2 = binfile[:, 2]
coor1_bins = np.linspace(min(coor1), max(coor1), num = Nsamples + 1)
coor2_bins = np.linspace(min(coor2), max(coor2), num = Nsamples + 1)
bins = np.linspace(0, len(coor1), num = Nsamples + 1)


#-------------Divide1--------------#
plt.figure()
plt.plot(coor1, coor2, 'k.')
Ncount = []
for i in range(Nsamples):
    plt.plot(coor[int(bins[i]):int(bins[i+1]), 0], coor[int(bins[i]):int(bins[i+1]), 1], '.', alpha = 0.5)
    Ncount.append(len(coor[int(bins[i]):int(bins[i+1]), 0]))
    np.savetxt('./Data_Output/miniJPAS_subregions/miniJPAS_' + bin_lum[0] + '_' + bin_lum[bin_num+1] + '_field_%i'%i + '_over%i'%Nsamples +'.txt', coor[int(bins[i]):int(bins[i+1]), :], fmt = "%.6f  %.6f",delimiter = "\n")
np.savetxt('./Data_Output/miniJPAS_subregions/miniJPAS_Ncount_' + bin_lum[0] + '_' + bin_lum[bin_num+1] + '_over%i'%Nsamples +'.txt', Ncount, fmt = "%.i")
plt.show();exit()


#-------------Divide2--------------#
plt.figure()
plt.plot(coor1, coor2, 'k.')
for i in range(Nsamples):
    coor1_remain = []; coor2_remain = []
    for j in range(len(coor1)):
        if coor1_bins[i] <= coor[j, 0] <= coor1_bins[i + 1] or coor2_bins[i] <= coor[j, 1] <= coor2_bins[i + 1]:
            coor1_remain.append(coor[j, 0])
            coor2_remain.append(coor[j, 1])
    plt.plot(coor1_remain, coor2_remain, '.', alpha = 0.5)
    np.savetxt('./Data_Output/miniJPAS_subregions/miniJPAS_' + bin_lum[0] + '_' + bin_lum[bin_num+1] + '_field_%i'%i + '_over%i'%Nsamples +'.txt', np.column_stack((coor1_remain, coor2_remain)), fmt = "%.6f  %.6f",delimiter = "\n")
plt.show()







