import numpy as np
from astropy.coordinates import *
import sys
import pandas as pd
from math import *

import colossus
from colossus import *
from colossus.cosmology import *
from colossus.halo import *
from colossus.lss import *

from scipy.interpolate import interp1d
import scipy.special as special
from scipy.integrate import quad

import Func_Cosmo as Cfunc
from param_fac import *

params = {'flat': True, 'H0': 100., 'Om0': 0.3089, 'Ob0': 0.0486, 'sigma8': 0.8161, 'ns': 0.9667}
cosmo = cosmology.setCosmology('miniJPAS_Cosmo', params)
## This code is to calculate Auto correlation power spectrum for miniJPAS.


#-------Global parameters---------#
global bin_lum
bin_lum = ['-25.2', '-17.5', '-20.0', '-21.0', '-22.0', '-23.0']


#-------------CODES---------------#
def get_z(ibin_num, zfilepath = './Data_Input/JPAS/', binfilepath = './Data_Output/'):
    zfile = (zfilepath + 'minijpas_pdr201912_cat_wAbsMag_v1.csv') 
    full_data = np.array(pd.read_csv(zfile))

    binfile = (binfilepath + 'miniJPAS_BB_LSS_iband_selection_' + bin_lum[ibin_num] + '.txt')
    bin_data = np.loadtxt(binfile)
    z_bin_data = np.zeros(len(bin_data))
    loss_id_z = []
    for i in range(len(z_bin_data)):
        if float(bin_lum[ibin_num]) <= bin_data[i, 16] <= float(bin_lum[ibin_num + 1]):
            try:
                id_z = (list(full_data[:, 0])).index(int(bin_data[i, 0]))
                z_bin_data[i] = full_data[id_z, 35]
            except:
                loss_id_z.append(int(bin_data[i, 0]))
                z_bin_data[i] = 100.
    np.savetxt(binfilepath + '/miniJPAS_' + bin_lum[ibin_num] + '_losszzIndex.dat', loss_id_z, '%i')
    np.savetxt(binfilepath + '/miniJPAS_' + bin_lum[ibin_num] + '_zz.dat', z_bin_data, '%4f')
    return z_bin_data


#---------------------------------#
def Winfunc_auto(izz, izz_dist = None, path = False, dndzopt = False):
## window function for galaxy-galaxy auto-correlation
## dndz here is from file.
    if path == None:
        raise Exception("Error: line", sys._getframe().f_lineno, ": Please input file path if you want dndz input from files, and file in a fixed format.")
    elif path == False:
        z_dist = izz_dist
    else:
        z_dist = np.load(path + 'z_ML.npy', encoding = 'latin1') 

    zbins = np.linspace(0.2, 0.6, num = 50)
    zbins_mid = np.zeros(len(zbins) - 1)
    dndz = np.zeros(len(zbins) - 1)
    dndz_deltazz = np.zeros(len(zbins) - 1)
    for i in range(len(zbins) - 1):
        zbins_mid[i] = (zbins[i] + zbins[i+1])/2.
        deltazz = (zbins[i+1] - zbins[i])
        dndz[i] = (((zbins[i]<z_dist) & (z_dist<zbins[i+1])).sum())
        dndz_deltazz[i] = dndz[i]/deltazz
    if dndzopt:
        return zbins_mid, dndz
    import scipy.integrate as integrate
    spl = interp1d(zbins_mid, dndz_deltazz/sum(dndz), kind = 'linear', fill_value = 'extrapolate')
    winres = spl(izz) 
    return winres


#---------------------------------#
def Shot_Noise_CAT():
    cnoise = 4. *np.pi* fsky / catntot
    return cnoise


#----------------------------#
def __dndm_func(imm, izz, iq_in = 'M', iq_out = 'dndlnM'):
    mfun = mass_function.massFunction(imm, izz, q_in = iq_in, q_out = iq_out, mdef = 'fof', model = 'watson13')/imm
    ## Unit: Mpc^-3 h^3
    return mfun


#----------------------------#
def N_hod(imm):
    N_cent = 0.5*(1. + special.erf( (np.log10(imm) - np.log10(iM_cutmin_ps)) /isigma_logm_ps))
    if imm > iM_0_ps:
        N_sat = N_cent*((imm - iM_0_ps) /iM_1_ps)**ialpha_hod_ps
    else:
        N_sat = 0.
    return N_cent, N_sat


#----------------------------#
def N_hod_3param(imm, iM_min_fix = False, iM_1_fix = False, ialpha_fix = False):
    if iM_min_fix:
        N_cent = np.where((imm-iM_min_fix) > 0., 1, 0)
        N_sat = N_cent*(np.maximum(imm - iM_min_fix, 0)/iM_1_fix)**ialpha_fix
    else:
        N_cent = np.where((imm-iM_min_ps) > 0., 1, 0)
        N_sat = (np.maximum(imm - iM_min_ps, 0)/iM_1_ps)**ialpha_ps
    return N_cent, N_sat


#----------------------------#
def Nbar(izz, iM_min, iM_max):
    mm_inte_simps = np.logspace(iM_min, iM_max, num = 30)
    N_cent, N_sat= N_hod_3param(mm_inte_simps)
    nbar_simps = __dndm_func(mm_inte_simps, izz) *(N_cent + N_sat)
    nbar = np.trapz(nbar_simps, mm_inte_simps, 30)
    ## Unit: Mpc^-3 h^3
    return nbar


#----------------------------#
def Ph_gg_forinte(imm, izz, ikk):
    ncent, nsat = N_hod_3param(imm)
    dndm = __dndm_func(imm, izz)
    ftrho = FtRho_NFW(imm, ikk, izz)
    bias = lss.bias.haloBias(imm, izz, mdef = '200c', model = 'tinker10')
    p1hgg = dndm *(2 *ncent*nsat*ftrho**2 + ncent*nsat**2 *ftrho**2)
    p2hgg = dndm *(ncent + nsat) *bias *ftrho
    return p1hgg, p2hgg


#----------------------------#
def __integrate_PS_gg(ikk, izz, iM_min, iM_1, ialpha):
    global iM_min_ps, iM_1_ps, ialpha_ps
    nn = 30
    iM_min_ps = iM_min *np.ones(nn); iM_1_ps = iM_1 *np.ones(nn); ialpha_ps = ialpha *np.ones(nn)

    mm_inte_simps = np.logspace(M_min_dm, M_max_dm, num = nn)
    ps_inte = 0.; ps_inte1 = 0.; ps_inte2 = 0.
    nbar2 = Nbar(izz, M_min_dm, M_max_dm)**2
    if nbar2 != 0:
        p1h_simps, p2h_simps = Ph_gg_forinte(mm_inte_simps, izz, ikk)
    ps_inte1 = np.trapz(p1h_simps, x = mm_inte_simps)/nbar2
    ps_inte2 = np.trapz(p2h_simps, x = mm_inte_simps)**2/nbar2
    ps_inte2 = ps_inte2 *Plin(ikk, izz)
    ps_inte = ps_inte1 + ps_inte2
    if verbose:
        print(' %.2e'%ikk, '%.4f'%izz, '%.3e'%(ps_inte1), '%.3e'%(ps_inte2), '%.3e'%(ps_inte))
    return ps_inte, ps_inte1, ps_inte2


#----------------------------#
def __integrate_CAPS(izz, ibin_num, iM_min, iM_1, ialpha):
    izz_dist = np.sort(np.loadtxt('./Data_Output/miniJPAS_redshift/miniJPAS_zz_nostar_iband_selection_' + bin_lum[0] + '_' + bin_lum[ibin_num+1] + '.txt')[:, 5])
    win_gg = Winfunc_auto(izz, izz_dist)
    ps = np.zeros(len(izz))
    for i in range(len(izz)):
        ps[i] = __integrate_PS_gg(ikk_inte[i], izz[i], iM_min, iM_1, ialpha)[0]
    caps_res = win_gg**2 *ps /Cfunc.Distfunc(izz, unit = 'Mpc')**2
    return caps_res


#----------------------------#
def __integration_Ndensity_g(imm, izz, iM_min, iM_1, ialpha):
    dncent, dnsat = N_hod_3param(imm, iM_min, iM_1, ialpha)
    nall = __dndm_func(imm, izz)*(dncent + dnsat)
    return nall


#----------------------------#
def __integration_Ndensity(imm, ibin_num, izz):
    Ndata = np.array([4887., 4257., 2817., 973., 130.])
    nall = __dndm_func(imm, izz)*Ndata[ibin_num + 1]
    return nall


#----------------------------#
def Ndensity_g(izz, iM_min, iM_1, ialpha):
    return quad(__integration_Ndensity_g, 1e10, 1e16, args=(izz, iM_min, iM_1, ialpha))[0]


#----------------------------#
def TheoryDATA(iM_min, iM_1, ialpha, izz_inte, iell, ibin_num):
    caps = np.zeros([len(iell)])
    for i in range(len(iell)):
        global ikk_inte
        ikk_inte = (iell[i]+0.5)/Cfunc.Distfunc(izz_inte, unit = 'Mpc')
        caps_inte = __integrate_CAPS(izz_inte, ibin_num, iM_min, iM_1, ialpha) \
                    /3e5 *Cfunc.Hzfun(izz_inte)
        caps[i] = np.trapz(caps_inte, x = izz_inte)
    return caps
    
    
#----------------------------#
def MeasurementDATA(ibin_num):
    ell_data = np.loadtxt('./Data_Output/miniJPAS_cl/miniJPAS_' + bin_lum[0] + '_' + bin_lum[ibin_num+1] + '_cl.dat')[2:, 0]
    Mdata = np.loadtxt('./Data_Output/miniJPAS_cl/miniJPAS_' + bin_lum[0] + '_' + bin_lum[ibin_num+1] + '_cl.dat')[2:, 1]
    Mvar = np.loadtxt('./Data_Output/miniJPAS_cl/miniJPAS_' + bin_lum[0] + '_' + bin_lum[ibin_num+1] +  '_error.dat')[2:]
    return ell_data, Mdata, Mvar


#----------------------------#
def MeasurementNdensity(ibin_num, izz):
    return quad(__integration_Ndensity, M_min, M_max, args=(ibin_num, izz))[0]


#----------------------------#
def Plin(ikk, izz):
## adiabatic matter power spectrum
    if ikk < 1.e20:
        plin = cosmo.matterPowerSpectrum(ikk, izz, model = 'eisenstein98')
    else:
        plin = 0.
    ## Unit: Mpc^3 h^-3
    return plin
    
    
#----------------------------#
def FtRho_NFW(imm, ikk, izz, forGas = False):
# Forier transformation for nfw rho, see <Func_Halo> for whole derivation.
    cc = 9. /(imm/1.e14*cosmo.h)**0.172/(1. + izz)
    #cc = concentration.concentration(imm/cosmo.h, '200c', izz, model = 'diemer19')
    if forGas:
        rvir = (imm/(4.*np.pi/3. *(cosmo.h**2 *200. *1.e9 *cosmo.rho_b(0.))))**(1./3)
    else:
        rvir = (imm/(4.*np.pi/3. *(cosmo.h**2 *200. *1.e9 *(cosmo.rho_m(0.) - cosmo.rho_b(0.)))))**(1./3)
    rs = rvir/cc

    si1=special.sici(ikk*rs)[0]
    ci1=special.sici(ikk*rs)[1]
    si2=special.sici(ikk*(rs+rvir))[0]
    ci2=special.sici(ikk*(rs+rvir))[1]
    si3=special.sici(ikk*rvir)[0]
    ci3=special.sici(ikk*rvir)[1]

    ftrho=rs**3*(np.sin(ikk*rs)*(si2-si1)-np.sin(ikk*rvir)/\
            (ikk*(rs+rvir))+np.cos(ikk*rs)*(ci2-ci1))\
            *(imm/(cosmo.rho_c(0.)*cosmo.h**2*1.e9)/rs**3/(np.log(1.+cc)-cc/(1.+cc)))
    ftrho = ftrho/imm*(cosmo.h**2*1.e9*cosmo.rho_c(0))
    ## Unit: Mpc^3 h^-3
    return ftrho



