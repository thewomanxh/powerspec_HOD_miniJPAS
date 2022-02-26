import numpy as np
import sys
from scipy.optimize import minimize, LinearConstraint
from scipy.interpolate import RegularGridInterpolator

## This code stores some functions to fit the HOD parameters.


#----------------------------#
def read_data(f = './Data_Output/Pablo_out/3param_out.csv', M1 = 0, M_min = 0, alpha = 0):
    data = np.loadtxt(f)
    if M1 == 1:
        logM1 = data[:, 6:10]
        return 10.**logM1
    if M_min == 1:
        logMmin = data[:, 2:6]
        return 10.**logMmin
    if alpha == 1:
        alpha_va = data[:, 10:14]
        return alpha_va


#----------------------------#
def read_Xiudata(f = './Data_Output/Xiu_out/hod_3param.dat', M1 = 0, M_min = 0, alpha = 0):
    data = np.loadtxt(f)
    if M1 == 1:
        logM1 = data[:, 4:7]
        #logM1[:, 1] += logM1[:, 0]; logM1[:, 2] += logM1[:, 0]
        return 10.**logM1
    if M_min == 1:
        logMmin = data[:, 1:4]
        #logMmin[:, 1] += logMmin[:, 0]; logMmin[:, 2] += logMmin[:, 0]
        return 10.**logMmin
    if alpha == 1:
        alpha_va = data[:, 7:10]
        return alpha_va


#----------------------------#
def alpha_model(magnitude, a, b, c):
    log_lum_z0 = 10**((4.76 - magnitude)/2.5)
    alpha_res = np.log10(c + (a*np.log10(log_lum_z0))**b)
    return alpha_res


#----------------------------#
def M0_model(magnitude, a, b):
    log_lum_z0 = 10**((4.76 - magnitude)/2.5)
    M0_res = 10.**(a*np.log10(log_lum_z0) + b)
    return M0_res
    

#----------------------------#
def M1_model_prior(logm, L_s, M_t, a_m):
    lum = L_s*(10**logm/M_t)**a_m * np.exp(1-(M_t/10**logm))
    magnitude = 4.76 - 2.5*np.log10(lum)
    return magnitude


#----------------------------#
def M1_model(magnitude, L_s, M_t, a_m):
    return 10.**M_interpolator(magnitude, L_s, M_t, a_m)


#----------------------------#
def Mmin_model():
    ## We do not have this parameter temporarily.
    return None

#----------------------------#
def sigma_model():
    ## We do not have this parameter temporarily.
    return None

#----------------------------#
def M_interpolator(magout, L_s, M_t, a_m):
    log_mass = np.arange(12, 16, 0.001)[::-1]
    magnitudes = M1_model_prior(log_mass, L_s, M_t, a_m)
    f = RegularGridInterpolator((magnitudes,), log_mass, bounds_error=False, fill_value=None)
    res = (f(magout))
    return res


#----------------------------#
def chi2(theta):
    if Pablo_fit == 1:
        ms_raw = read_data(f = './Data_Output/Pablo_out/3param_out.csv', M1 = M1_fit, M_min = M_min_fit, alpha = alpha_fit)
        mag_raw = np.loadtxt('./Data_Output/Pablo_out/3param_out.csv')[:-1, 0]
        ms_data = ms_raw[:-1, 2]
        ms_var = ms_raw[:-1, 3] - ms_raw[:-1, 1]
    if Xiu_fit == 1:
        ms_raw = read_Xiudata(f = './Data_Output/Xiu_out/hod_3param.dat', M1 = M1_fit, M_min = M_min_fit, alpha = alpha_fit)
        mag_raw = np.loadtxt('./Data_Output/Xiu_out/hod_3param.dat')[:, 0]
        ms_data = ms_raw[:, 0]
        ms_var = ms_raw[:, 2] - ms_raw[:, 1]
    
    if M_min_fit == 1:
        th_data = M0_model(mag_raw, theta[0], theta[1])
    if M1_fit == 1:
        th_data = M1_model(mag_raw, theta[0], theta[1], theta[2])
    if alpha_fit == 1:
        th_data = alpha_model(mag_raw, theta[0], theta[1], theta[2])
    Fitsdata = th_data - ms_data
    chi2 = sum(Fitsdata**2 / ((ms_var))**2)
    return chi2


#----------------------------#
def Fit_func(M1 = 0, M_min = 0, alpha = 0, Pablo = 0, Xiu = 0):
    global M1_fit, M_min_fit, alpha_fit, Pablo_fit, Xiu_fit
    M1_fit = M1; M_min_fit = M_min; alpha_fit = alpha; Pablo_fit = Pablo; Xiu_fit = Xiu
    if M_min == 1:
        pos = [1.78353177, -5.98024172]
        bounds_prior = ((0., 3.), (-10., 0.))
        imethod = 'Nelder-Mead'
    if M1 == 1:
        pos = [3.70558341e+09, 4.77678420e+11, 3.05963831e-01]
        bounds_prior = ((1.e8, 1.e10), (1.e11, 1.e14), (0.1, 3.))
        imethod = 'TNC'
    if alpha == 1:
        pos = [0.0982841, 80.27598, 10.0]
        bounds_prior = ((0., 1.), (70., 150.), (1., 50.))
        imethod = 'TNC'
    soln = minimize(chi2, pos, method = imethod, bounds = bounds_prior, tol = 1.e-7, options={'disp': True, 'maxiter': 1000})
    return soln.x


