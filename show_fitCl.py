import numpy as np
import matplotlib.pyplot as plt

from plotParam import *
from param_fac import *
import powerspec_HOD as HOD
import time

## This code is to show the best fit parameters compared with the data.


#-----Global parameters------#
global bin_lum, ibin_num
bin_lum = ['-25.2', '-17.5', '-20.0', '-21.0', '-22.0', '-23.0']
ibin_num = 2

global M_1_fit, M_min_fit, alpha_fit
#zz_inte = np.linspace(0.001, 1.5, num = 20, endpoint = True)
zz_inte = np.linspace(0.2, 0.6, num = 50, endpoint = True)

global Ndata, Nerr
Ndata = np.array([0.013538, 0.01167, 0.008965, 0.005926, 0.003093, 0.001254])
Nerr = np.array([0.000594, 0.000521, 0.000454, 0.000315, 0.000218, 0.000113])

TheoryOnly = False
t1 = time.time()

#----------------------------#
ell_data, ms_data, ms_var = HOD.MeasurementDATA(ibin_num)
ell_data = ell_data[2:48]
ms_data = ms_data[2:48]
ms_var = ms_var[2:48]
Fitting = np.array([[11.501961, 13.198678, 1.089903], [11.573894, 13.240629, 1.130207], [11.698694, 13.318858, 1.135367], [11.882649, 13.522107, 1.135342], [12.167323, 13.986156, 0.827222]])
#Fitting = np.array([[11.454, 13.10249, 0.534], [11.660, 12.916, 0.715], [11.711, 13.152, 0.789], [12.083, 13.547, 0.899], [12.146, 13.986, 0.701]])
#Fitting = np.array([[11.245, 11.224, 0.500], [11.661, 13.007, 0.655], [11.865, 12.838, 0.503], [11.995, 12.929, 0.554], [12.146, 13.986, 0.701]])
## chi2:  1508, 71, 399, 75 

M_min_fit = 10.**(Fitting[ibin_num, 0])
M_1_fit = 10**(Fitting[ibin_num, 1])
alpha_fit = Fitting[ibin_num, 2]


if TheoryOnly == False:
    print('---Best Fits for five parameters: %.4e'%M_min_fit, ' %.4e'%M_1_fit, ' %.4f'%alpha_fit, '\n', ' ')
    #print('---Measurement Data is: \n', ms_data, '\n', ' ')

    Tout = HOD.TheoryDATA(M_min_fit, M_1_fit, alpha_fit, zz_inte, ell_data, ibin_num)
    Nall = HOD.Ndensity_g(0.4, M_min_fit, M_1_fit, alpha_fit)
    Nfit = Nall - Ndata[ibin_num]
    #print('---Theoretical Data is: \n', Tout)
    print('\n---Total numbers: ', Nall, Ndata[ibin_num])
    print('\n---chi2: ', sum((ms_data-Tout)**2/ms_var**2) + Nfit**2/Nerr[ibin_num]**2)
    print('\n---time: ', time.time()-t1)


#----------------------------#
if TheoryOnly:
    colors = ['#b7e1a1', '#75b84f', '#699d4c', '#3f9b0b', '#388004']
    fig = plt.figure(1, (18, 14))
    plt.loglog()

    ell_data, ms_data, ms_var = HOD.MeasurementDATA(0)
    for i in range(len(bin_lum) - 1):
        Tout = HOD.TheoryDATA(10.**(Fitting[i, 0]), 10**(Fitting[i, 1]), Fitting[i, 2], zz_inte, ell_data, i)
        plt.plot(ell_data, Tout, color = colors[i], label = bin_lum[i+1])

    plt.xlim(min(ell_data)-10, max(ell_data)+1000)
    plt.ylim(1.e-9, 1.e-5)
    plt.xlabel(r'$\mathrm{Multipole\, }\ell$', fontsize = 14)
    plt.ylabel(r'$C_\ell$', fontsize = 15)
    plt.legend(loc = 'best', fontsize = 12, frameon = False)
    plt.tight_layout()
    plt.savefig('./Data_Output/Fit_res/All_Cl_miniJPAS.png')
    plt.show()

else:
    fig = plt.figure(1, (18, 14))
    plt.loglog()
    plt.errorbar(ell_data, ms_data, yerr = ms_var, color = 'blue', fmt = 'h', linewidth = 1, markersize = 7, capsize = 5, label = 'Measurment')
    plt.plot(ell_data, Tout, label = 'Theoretical fit')
    plt.xlim(min(ell_data)-10, max(ell_data)+1000)
    plt.ylim(1.e-9, 1.e-5)
    plt.text(min(ell_data)+20, 2.e-8, s = str(Fitting[ibin_num, :]), fontsize = 12)
    plt.xlabel(r'$\mathrm{Multipole\, }\ell$', fontsize = 14)
    plt.ylabel(r'$C_\ell$', fontsize = 15)
    plt.legend(loc = 'best', fontsize = 14, frameon = False)
    plt.tight_layout()
    plt.savefig('./Data_Output/Fit_res/CompareCl_miniJPAS_' + bin_lum[0] + '_' + bin_lum[ibin_num+1] + '.png')
    plt.show()


