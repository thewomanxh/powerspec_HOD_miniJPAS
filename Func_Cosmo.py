import numpy as np
import os
import sys
from math import *
import scipy.integrate as integrate
from scipy.interpolate import interp1d, interp2d

import colossus
from colossus import *
from colossus.cosmology import *

#cosmo = cosmology.setCosmology('planck18')
params = {'flat': True, 'H0': 67.74, 'Om0': 0.3089, 'Ob0': 0.0486, 'sigma8': 0.8161, 'ns': 0.9667}#, 'h': 0.6774}
cosmo = cosmology.setCosmology('miniJPAS_Cosmo', params)

## This is the code for cosmological quantities.


#--------------CODES---------------#
def Distfunc(izz, zi = 1.e-50, form = 'comoving', unit = 'cm', ih = False):
## This is the cosmological distance after integate, choose distance mode and units.
    if form == 'luminosity':
        dist = cosmo.luminosityDistance(izz)
    elif form == 'comoving':
        dist = cosmo.comovingDistance(zi, izz)
    elif form == 'angularDiameter':
        dist = cosmo.angularDiameterDistance(izz)
    ## all in units of Mpc/h
    else:
        raise Exception("Error: line", sys._getframe().f_lineno, ": Please input recognized format.")
    if ih == False:
        ih = cosmo.h
    if unit == 'cm':
        res = dist *3.08567758128E+24 *ih
    elif unit == 'Mpc':
        res = dist *ih
    else:
        raise Exception("Error: line", sys._getframe().f_lineno, ": Please input recognized format.")
    return res


#----------------------------------#
def Dist_comov_inte(izz):
    d = 3.e5/Hzfun(izz)
    ## Unit: Mpc, 3e5 is light speed in km/s.
    return d


#----------------------------------#
def Ezfunc(izz):
## the same with hzfun in angular.
    f = np.sqrt(cosmo.Ode(0) + cosmo.Om(0) *(1. + izz)**3 + \
            (1 - cosmo.Om(0) - cosmo.Ode(0)) *(1. + izz)**2)
    return f
 
 
#---------------------------------#
def Hzfun(izz):
## H0 unit is km/s/Mpc.
    Hz = cosmo.H0 *Ezfunc(izz)
    return Hz


#---------------------------------#
def Timefunc(izz, form = 'UniAge', unit = 's', iinverse = False):
## if inverse for look back time is true, input should be time (in Gyr), and output is redshift. NOTE: now the maximum izz is <500.
    if form == "UniAge":
        tt = cosmo.age(izz)
    elif form == "lookbacktime":
        tt = cosmo.lookbackTime(izz, inverse = iinverse)
    elif form == "hubbletime":
        tt = cosmo.hubbleTime(izz)
    ## all units in Giga years
    else:
        raise Exception("Error: line", sys._getframe().f_lineno, ": Please input recognized format.")
    if unit == 's':
        res = tt *1.e9 *365 *24 *3600.
    elif unit == 'year':
        res = tt *1.e9
    elif unit == 'Gyr':
        res = tt
    else:
        raise Exception("Error: line", sys._getframe().f_lineno, ": Please input recognized format.")

    return res


#----------------------------------#
def Time2zz(itt, unit = 's'):
## This function is from universe age to redshift, redshift range is 1.e-4 to 1.e6.
    ntt = 5000
    zzin = np.linspace(0., 1.e6, num = ntt, endpoint = True)

    if os.path.exists('./Data_Output/_table_out/zz_age_interpolate.dat') == False:
        print("---Calculating redshift VS age, total ", ntt, " templates:")
        ttout = np.zeros(ntt)
        for i in range(ntt):
            sys.stdout.write('  {0}\r'.format(i+1))
            sys.stdout.flush()
            ttout[i] = zz2Time(zzin[i])
        np.savetxt('./Data_Output/_table_out/zz_age_interpolate.dat', np.column_stack([ttout, zzin]), fmt = '%.8e    %.8e')
    else:
        filepf = np.loadtxt('./Data_Output/_table_out/zz_age_interpolate.dat')
        ttout = filepf[:, 0]
        zzin = filepf[:, 1]
    func = interp1d(ttout, zzin, kind = 'linear')

    if unit == 's':
        tt = itt
    elif unit == 'year':
        tt = itt *365 *24 *3600
    elif unit == 'Gyr':
        tt = itt *1.e9 *365 *24 *3600
    else:
        raise Exception("Error: line", sys._getframe().f_lineno, ": Please input recognized format.")
    tt = max(tt, ttout[-1])
    zzres = func(tt)

    return zzres
 

#----------------------------------#
def zz2Time(izz, unit = 's'):
## This is calculating the universe age of the input redshift izz, the redshift range can be any number.
    #tt = integrate.quad(inte_time, 0., 1100.)[0] - integrate.quad(inte_time, 0., izz)[0]
    time_inte = lambda a : 1./Hzfun(a)/(1.+a)/3.240779e-20
    ## Unit: s, 1/3.240779e-20. is Mpc to km.
    tt = integrate.quad(time_inte, 0., np.inf)[0] - integrate.quad(time_inte, 0., izz)[0]

    if unit == 's':
        return tt
    elif unit == 'year':
        return tt /(365. *24. *3600.)
    elif unit == 'Gyr':
        return tt /(1.e9 *365 *24 *3600)
    else:
        raise Exception("Error: line", sys._getframe().f_lineno, ": Please input recognized format.")


#----------------------------------#
def dVdz_inte(izz):
## This is the derivation comoving volume for computing window
    res = Dist_comov_inte(izz)**3 *4*np.pi/3
    ## Unit: Mpc^3 sr^-1
    return res


#----------------------------------#
def Vz(izz, unit = 'cm'):
## This is the Volume comoving.
    res = integrate.quad(dVdz_inte, izz, np.inf)[0] *4.*np.pi
    if unit == 'cm':
        res = res *pow((3.08567758128e24), 3)
    elif unit == 'Mpc':
        pass
    else:
        raise Exception("Error: line", sys._getframe().f_lineno, ": Please input recognized format.")
    ## Unit: Mpc^3
    return res


#---------------------------------#
def dndefunc_DM(iee):
    """
    It's from PPPC table to calculate dark matter dnde spectrum.
    Choosing channel by the last item.
    mass for annihilation, if you calculate decay mode, mass should be the half.
    The process of dark matter anni or decay is included EW and at production gamma.
    Pay attention with the mass and energy unit, GeV.
    
    """
    dndedata = np.loadtxt('./Data_Input/_table_in/AtProduction_gammas.dat')
    listname=['mdm', 'log10x', 'eL', 'eR', 'e', 'mul', 'mur', 'mu', 'taol', 'taor',\
            'tao', 'q', 'c', 'b', 't', 'WL', 'WT', 'W', 'ZL', 'ZT', 'Z', 'g', 'gamma',\
            'h', 'nue', 'numu', 'nutao', 'v_e', 'v_nu', 'v_tao']
    idx = listname.index(channel)
    log10x = []
    dnde = []
    #mDMdec = mDMdec/2.
    for i in range(2, len(dndedata)):
        if dndedata[i, 0] == mDMdec/2.:
            log10x.append(dndedata[i, 1])
            dnde.append(dndedata[i, idx])
    log10x = np.array(log10x)
    dnde = np.array(dnde)
    EE = np.array([np.log10(iee /mDMdec/2.)])
    for i in range(len(EE)):
        func1 = interp1d(log10x, dnde, kind = 'quadratic')
        func2 = extrap1d(func1)
        dnde_new = func2(EE)
    return dnde_new


#---------------------------------#
def Optfunc(izz, iee):
## Optical data interpolate and extrapolate, input iee in GeV unit.
    optdata = np.loadtxt('./Data_Input/_table_in/opdep_fixed.txt') ##GeV
    z_opt=np.array([0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,\
            0.09,0.10,0.11,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,\
            0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.2,1.4,\
            1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,\
            4.2,4.4,4.6,4.8,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0])
    ## because we just need low-z, so throw z>6 away.
    z_opt = z_opt[:51]
    e_opt = optdata[:,0]/1000.
    opt = optdata[:,1:52]
    newfunc = interp2d(z_opt,e_opt,opt,kind='cubic')

    if isinstance(iee, float) == False and isinstance(izz, float) == False:
        opt_new=np.zeros((len(iee), len(izz)))
        ## EE and zz is list.
        for i in range(len(iee)):
            for j in range(len(izz)):
                for k in range(len(z_opt)):
                    for l in range(len(e_opt)):
                        if izz[j] <= z_opt[k] and iee[i] <= e_opt[l] and opt[l,k] <= 0.:
                            opt_new[i,:j] = 0.
                            continue
                        else:
                            opt_new[i,j] = newfunc(izz[j], iee[i])
                            if opt_new[i,j] < 0.:
                                opt_new[i,j] = 0.

    elif isinstance(izz,float) == False and isinstance(iee,float) == True:
        opt_new = np.zeros((len(izz)))
        for j in range(len(izz)):
            for k in range(len(z_opt)):
                for l in range(len(e_opt)):
                    if izz[j] <= z_opt[k] and iee <= e_opt[l] and opt[l,k] <= 0.:
                        opt_new[j] = 0.
                        continue
                    else:
                        opt_new[j] = newfunc(izz[j], iee)
                        if opt_new[j] < 0.:
                            opt_new[j] = 0.

    else:
        for k in range(len(z_opt)):
            for l in range(len(e_opt)):
                if izz <= z_opt[k] and iee <= e_opt[l] and opt[l,k] <= 0.:
                    opt_new = 0.
                else:
                    opt_new = newfunc(izz, iee)
                    if opt_new < 0.:
                        opt_new = 0.
    return opt_new


