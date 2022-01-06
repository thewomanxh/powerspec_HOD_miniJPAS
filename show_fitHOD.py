import numpy as np
import fit_HOD as hodfit
import matplotlib.pyplot as plt
import sys
from plotParam import * 

## This file is to show the fit HOD parameters results, this results are supposed to hodpy_miniJPAS/lookup/hod_params.dat.


M_min_set = 1
M1_set = 0
alpha_set = 0

Pablo_set = 1
Xiu_set = 0


#----------------------------#
mag = np.arange(-23.5, -17., 0.1)[::-1]
if Pablo_set == 1:
    ms_raw = hodfit.read_data(f = './Data_Output/Pablo_out/3param_out.csv', M1 = M1_set, M_min = M_min_set, alpha = alpha_set)
    mag_raw = np.loadtxt('./Data_Output/Pablo_out/3param_out.csv')[:-1, 0]
    save_name = '_fit_Pablo.png'
    ms_data = ms_raw[:-1, 2]
    
if Xiu_set == 1:
    ms_raw = hodfit.read_Xiudata(f = './Data_Output/Xiu_out/hod_3param.dat', M1 = M1_set, M_min = M_min_set, alpha = alpha_set)
    mag_raw = np.loadtxt('./Data_Output/Xiu_out/hod_3param.dat')[:-1, 0]
    save_name = '_fit_xiu.png'

M_3params = hodfit.Fit_func(M_min = M_min_set, M1 = M1_set, alpha = alpha_set, Pablo = Pablo_set, Xiu = Xiu_set)

if M_min_set == 1:
    Mmin_th = hodfit.M0_model(mag, M_3params[0], M_3params[1])
if M1_set == 1:
    M1_th = hodfit.M1_model(mag, M_3params[0], M_3params[1], M_3params[2])
if alpha_set == 1:
    alpha_th = hodfit.alpha_model(mag, M_3params[0], M_3params[1], M_3params[2])

if Pablo_set == 1:
    ms_data = ms_raw[:-1, 2]
    ms_var = ms_raw[:-1, 3] - ms_raw[:-1, 1]
if Xiu_set == 1:
    ms_data = ms_raw[:-1, 0]
    ms_var = ms_raw[:-1, 2] - ms_raw[:-1, 1]
    

#----------------------------#
fig = plt.figure(1, (18, 16))
plt.xlim([-17., -23.5])
plt.xlabel(r'$^{0.4}M_i - 5\log{h}$')
#plt.errorbar(mag_raw, ms_data, yerr = ms_var/2. + ms_raw[:-1, 0], fmt = '.', elinewidth = 1, capsize = 2, capthick = 0.5)
plt.errorbar(mag_raw, ms_data, yerr = ms_var, fmt = '.', elinewidth = 1, capsize = 2, capthick = 0.5)
if M_min_set == 1:
    plt.semilogy()
    plt.plot(mag, Mmin_th)
    plt.ylim(1.e10, 1.e13)
    plt.ylabel(r'$M_{min}(h^{-1}M_\odot)$')
    plt.savefig('./Data_Output/HOD_params_res/Mmin' + save_name)
if M1_set == 1:
    plt.semilogy()
    plt.plot(mag, M1_th)
    plt.ylim([1.e12, 1.e16])
    plt.ylabel(r'$M_1(h^{-1}M_\odot)$')
    plt.savefig('./Data_Output/HOD_params_res/M1' + save_name)
if alpha_set == 1:
    plt.plot(mag, alpha_th)
    plt.ylim([0.2, 1.6])
    plt.ylabel(r'$\alpha$')
    plt.savefig('./Data_Output/HOD_params_res/alpha' + save_name)
#plt.tight_layout()
plt.show()

