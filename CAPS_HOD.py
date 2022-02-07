import os
import sys
import time
from scipy.optimize import minimize
import numpy as np
from scipy.integrate import quad
from colossus import *
from colossus.cosmology import *
from colossus.lss import *
import compute_HOD as HOD
import Func_Cosmo as Cfunc
from param_fac import *
params = {'flat': True, 'H0': 100., 'Om0': 0.3089, 'Ob0': 0.0486, 'sigma8': 0.8161, 'ns': 0.9667, 'h': 0.6774}
cosmo = cosmology.setCosmology('myCosmo', params)



#-----Global parameters------#
global bin_lum
bin_lum = ['-25.2', '-17.5', '-20.0', '-21.0', '-22.0', '-23.0']


#----------------------------#
def MeasurementDATA(ibin_num):
    ell_data = np.loadtxt('./Data_Output/miniJPAS_cl/miniJPAS_' + bin_lum[0] + '_' + bin_lum[ibin_num+1] + '_cl.dat')[2:, 0]
    Mdata = np.loadtxt('./Data_Output/miniJPAS_cl/miniJPAS_' + bin_lum[0] + '_' + bin_lum[ibin_num+1] + '_cl.dat')[2:, 1]
    Mvar = np.loadtxt('./Data_Output/miniJPAS_cl/miniJPAS_' + bin_lum[0] + '_' + bin_lum[ibin_num+1] +  '_error.dat')[2:]
    return ell_data, Mdata, Mvar


#----------------------------#
def TheoryDATA(iM_min, iM_1, ialpha, izz, iell, ibin_num):
    izz_dist = np.sort(np.loadtxt('./Data_Output/miniJPAS_redshift/miniJPAS_zz_nostar_iband_selection_' + bin_lum[0] + '_' + bin_lum[ibin_num+1] + '.txt')[:, 5])
    win_gg = HOD.Winfunc_auto(izz, izz_dist)

    ps = np.zeros([len(iell), len(izz)])
    caps_res = np.zeros([len(iell), len(izz)])
    caps = np.zeros([len(iell)])
    Comov_dist = Cfunc.Distfunc(izz, unit = 'Mpc')

    for j in range(len(iell)):
        for i in range(len(izz)):
            ps[j, i] = HOD.Integrate_PS_gg((iell[j]+0.5)/Comov_dist[i], izz[i], iM_min, iM_1, ialpha)[0]
        caps[j] = np.trapz( pow(win_gg, 2) *ps[j, :] /pow(Comov_dist, 2) /3e5 *Cfunc.Hzfun(izz), x = izz)
    return caps


#----------------------------#
def __integration_Ndensity_g(imm, izz, iM_min, iM_1, ialpha):
    dncent, dnsat = HOD.N_hod_3param(imm, iM_min, iM_1, ialpha)
    nall = (dncent + dnsat) *HOD.dndm_func(imm, izz)
    return nall


#----------------------------#
def Ndensity_g(izz, iM_min, iM_1, ialpha):
    imm_inte =  np.logspace(10., 16., num = 30)
    inte_Ndensity = __integration_Ndensity_g(imm_inte, izz, iM_min, iM_1, ialpha)
    res = np.trapz(inte_Ndensity, imm_inte)
    return res


#----------------------------#
def Volume_field(area, izzmin, izzmax):
    ## area given in steradians.
    dist_min = Cfunc.Distfunc(izzmin, unit = 'Mpc')
    dist_max = Cfunc.Distfunc(izzmax, unit = 'Mpc')
    volume_res = (area/3.) *(pow(dist_max, 3) - pow(dist_min, 3))
    return volume_res


#----------------------------#
def __integration_ms_Ndensity(imm, ibin_num, izz):
    Ndata = np.array([4887., 4257., 2817., 973., 130.])
    nall = HOD.dndm_func(imm, izz)*Ndata[ibin_num+1]
    return nall


#----------------------------#
def ms_Ndensity(ibin_num, izz):
    return quad(__integration_ms_Ndensity, 1e10, 1e16, args=(ibin_num, izz))[0]


