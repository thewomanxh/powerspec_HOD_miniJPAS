import numpy as np
import matplotlib.pyplot as plt
import fit_HOD as hodfit
from plotParam import *
import sys

if sys.argv[1] == 'M1':
    iM1 = 1; iM_min = 0; ialpha = 0
    ylims = [1.e11, 1.e15]
    yname = r'$M_{1}(h^{-1}M_\odot)$'
if sys.argv[1] == 'Mmin':
    iM_min = 1; iM1 = 0; ialpha = 0
    ylims = [1.e11, 1.e13]
    yname = r'$M_{min}(h^{-1}M_\odot)$'
if sys.argv[1] == 'alpha':
    ialpha = 1; iM1 = 0; iM_min = 0
    ylims = [0.5, 2.]
    yname = r'$\alpha$'


ms_pab = hodfit.read_data(f = './Data_Output/Pablo_out/3param_out.csv', M1 = iM1, M_min = iM_min, alpha = ialpha)
mag_pab = np.loadtxt('./Data_Output/Pablo_out/3param_out.csv')[:, 0]
pab_data = ms_pab[:, 2]
pab_var = ms_pab[:, 3] - ms_pab[:, 1]

ms_xiu = hodfit.read_Xiudata(f = './Data_Output/Xiu_out/hod_3param_output.dat', M1 = iM1, M_min = iM_min, alpha = ialpha)
mag_xiu = np.loadtxt('./Data_Output/Xiu_out/hod_3param_output.dat')[:, 0]
xiu_data = ms_xiu[:, 0]
if iM1 == 1 or iM_min == 1:
    xiu_var = xiu_data*(ms_xiu[:, 2]/np.log10(xiu_data)) - xiu_data*(ms_xiu[:, 1]/np.log10(xiu_data))
else:
    xiu_var = 10**(np.log10(xiu_data)+ms_xiu[:, 2]) - 10**(np.log10(xiu_data)-ms_xiu[:, 1])

#----------------------------#
fig = plt.figure(1, (18, 16))
plt.semilogy()
plt.xlim([-19., -23.])
plt.xlabel(r'$^{0.4}M_i - 5\log{h}$')
plt.errorbar(mag_pab, pab_data, yerr = pab_var, fmt = '.', elinewidth = 1, capsize = 2, capthick = 0.5)
plt.errorbar(mag_xiu, xiu_data, yerr = xiu_var, fmt = '.', elinewidth = 1, capsize = 2, capthick = 0.5)
plt.ylim(ylims)
plt.ylabel(yname)
plt.tight_layout()
plt.savefig('./Data_Output/Fit_res/param_comparison_' + str(sys.argv[1]) + '.png')
plt.show()


