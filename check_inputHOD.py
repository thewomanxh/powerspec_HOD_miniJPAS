import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator


Mmin_plot = False
M1_plot = False
alpha_plot = True

#------------Functions------------#
def lum2mag(luminosity):
    """
    Convert luminosity to absolute magnitude
    Args:
        luminosity: array of luminsities [Lsun/h^2]
    Returns:
        array of absolute magnitude [M-5logh]
    """
    return 4.76 - 2.5*np.log10(luminosity)

#---------------------------------#
def M1_model_prior(logm, L_s, M_t, a_m):
    lum = L_s*(10**logm/M_t)**a_m * np.exp(1-(M_t/10**logm))
    magnitude = 4.76 - 2.5*np.log10(lum)
    return magnitude
    
#---------------------------------#
def M_interpolator(magout, L_s, M_t, a_m):
    log_mass = np.arange(12, 16, 0.001)[::-1]
    magnitudes = M1_model_prior(log_mass, L_s, M_t, a_m)
    f = RegularGridInterpolator((magnitudes,), log_mass, bounds_error=False, fill_value=None)
    res = (f(magout))
    return res

#---------------------------------#
def M_function(log_mass, L_s, M_t, a_m):
    """
    Function describing how HOD parameters Mmin and M1 vary with magnitude.
    This function returns magnitude as a function of mass, and has the same 
    functional form as Eq. 11 from Zehavi 2011.

    Args:
        log_mass: array of the log10 of halo mass (Msun/h)
        L_s:      median luminosity of central galaxies (Lsun/h^2)
        M_t:      transition halo mass (Msun/h)
        a_m:      high-mass power law index     
    Returns:
        array of absolute magnitude
    """
    lum = L_s*(10**log_mass/M_t)**a_m * np.exp(1-(M_t/10**log_mass))
    magnitudes = lum2mag(lum)
    return magnitudes

#---------------------------------#
def mag2lum(magnitude):
    """
    Convert absolute magnitude to luminosity
    Args:
        magnitude: array of absolute magnitudes [M-5logh]
    Returns:
        array of luminosity [Lsun/h^2]
    """
    return 10**((4.76 - magnitude)/2.5)

#---------------------------------#
def alpha_func(abs_mag, alpha_A, alpha_B, alpha_C):
    """
    From line 435 in hod_bgs.py
    """
    
    log_lum = np.log10(mag2lum(abs_mag))
    
    return np.log10(alpha_C + (alpha_A*log_lum)**alpha_B)

#---------------------------------#
def Mmin(magnitude, redshift, f=1.0):
    """
    HOD parameter Mmin, which is the mass at which a halo has a 50%
    change of containing a central galaxy brighter than the magnitude 
    threshold
    Args:
        magnitude: array of absolute magnitude threshold
        redshift:  array of halo redshifts
    Returns:
        array of Mmin
    """
    return 10**M_interpolator(magnitude)

#---------------------------------#
def M1(magnitude, L_s, M_t, a_m, f=1.0):
    """
    HOD parameter M1, which is the mass at which a halo contains an
    average of 1 satellite brighter than the magnitude threshold
    Args:
        magnitude: array of absolute magnitude threshold
        redshift:  array of halo redshifts
    Returns:
        array of M1
    """
    return 10**M_interpolator(magnitude, L_s, M_t, a_m)

#---------------------------------#
def M0(magnitude, Mmin_A, Mmin_B):
    """
    HOD parameter M0, which sets the cut-off mass scale for satellites
    satellites
    Args:
        magnitude: array of absolute magnitude threshold
        redshift:  array of halo redshifts
    Returns:
        array of Mmin, since actually our M_min is M0
    """
    log_lum_z0 = (4.76 - magnitude)/2.5

    return 10**(Mmin_A*log_lum_z0 + Mmin_B)


#---------------------------------#
## Read in hod_params.dat
hod_params_np = np.loadtxt("../hodpy_miniJPAS/lookup/hod_params_pablo.dat", delimiter=',')

## Read in miniJPAS results
df_miniJPAS_pablo = np.loadtxt("./Data_Output/Pablo_out/3param_out.dat")

Mmin_A, Mmin_B = hod_params_np[2,:2]
M1_Ls, M1_Mt, M1_am = hod_params_np[1,:3]
alpha_A, alpha_B, alpha_C = hod_params_np[3,:3]

absmag_array = np.linspace(-19, -23, 100)
Mmin_array = M0(absmag_array, Mmin_A, Mmin_B)
M1_array = M1(absmag_array, M1_Ls, M1_Mt, M1_am)
alpha_array = alpha_func(absmag_array, alpha_A,alpha_B,alpha_C)


#--------------Plot---------------#
if alpha_plot:
    f, ax = plt.subplots()
    ax.plot(df_miniJPAS_pablo[:, 0], df_miniJPAS_pablo[:, 10], color='C0', label=r'HOD model - $w_p(r_p)$')
    ax.fill_between(df_miniJPAS_pablo[:, 0], df_miniJPAS_pablo[:, 11],df_miniJPAS_pablo[:, 13],color='C0', alpha=0.5)
    ax.set_xlim(-19.8, -22.2)
    ax.plot(absmag_array, alpha_array, color='C2', label='hodpy input')
    ax.set_ylim(0.4,1.4)
    ax.legend(loc='best')
    ax.set_xlabel(r'$M_r$ (limit)')
    ax.set_ylabel(r'$\alpha$')
    plt.show()

if Mmin_plot:
    f, ax = plt.subplots()
    ax.set_xlim(-19.8, -22.2)
    ax.plot(df_miniJPAS_pablo[:, 0], df_miniJPAS_pablo[:, 2], color='C0', label=r'HOD model - $w_p(r_p)$')
    ax.fill_between(df_miniJPAS_pablo[:, 0], df_miniJPAS_pablo[:, 3], df_miniJPAS_pablo[:,5],color='C0', alpha=0.5)
    ax.plot(absmag_array, np.log10(Mmin_array), color='C2', label='hodpy input')
    ax.legend(loc='best')
    ax.set_xlabel(r'$M_r$ (limit)')
    ax.set_ylabel(r'$\log_{10}(M_{\rm min})$')
    plt.show()

if M1_plot:
    f, ax = plt.subplots()
    ax.set_xlim(-19.8, -22.2)
    ax.plot(df_miniJPAS_pablo[:, 0], df_miniJPAS_pablo[:, 6], color='C0', label=r'HOD model - $w_p(r_p)$')
    ax.fill_between(df_miniJPAS_pablo[:, 0], df_miniJPAS_pablo[:, 7], df_miniJPAS_pablo[:,9],color='C0', alpha=0.5)
    ax.plot(absmag_array, np.log10(M1_array), color='C2', label='hodpy input')
    ax.set_xlabel(r'$M_r$ (limit)')
    ax.legend(loc='best')
    ax.set_ylabel(r'$\log_{10}(M_{\rm 1})$')
    plt.show()
