import os
import numpy as np
import time
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from multiprocessing import Pool
import corner
os.environ["OMP_NUM_THREADS"] = "1" 

import powerspec_HOD as HOD
from param_fac import *


#-----Global parameters------#
global bin_lum, ibin_num
bin_lum = ['-25.2', '-20.0', '-20.5', '-21.0', '-21.5', '-22.0', '-22.5', '-23.0']
#bin_lum = ['-25.2', '-17.5', '-20.0', '-21.0', '-22.0', '-23.0']
ibin_num = 4#int(sys.argv[1])


global M_1_fit, M_min_fit, alpha_fit
zz_inte = np.linspace(0.2, 0.6, num = 20, endpoint = True)

#global Ndata, Nerr
Ndata = np.array([0.013538, 0.01167, 0.008965, 0.005926, 0.003093, 0.001254])
Nerr = np.array([0.000594, 0.000521, 0.000454, 0.000315, 0.000218, 0.000113])

global ell_data, ms_data, ms_var, ell_interp
ell_data, ms_data, ms_var = HOD.MeasurementDATA(ibin_num)
ell_data = ell_data[0:93]
ms_data = ms_data[0:93]
ms_var = ms_var[0:93]
ell_interp = np.logspace(np.log10(200), np.log10(20000), 10) 
## FIXME: notice the range of \ell.


#----------------------------#
def log_likelihood(theta):
    M_min_fit, M_1_fit, alpha_fit = theta
    th_interp = HOD.TheoryDATA(10.**M_min_fit, 10.**M_1_fit, alpha_fit, zz_inte, ell_interp, ibin_num)
    func = interp1d(ell_interp, th_interp, kind = 'linear')
    th_data = func(ell_data)
    Fitsdata = th_data - ms_data

    #Nall = HOD.Ndensity_g(0.4, 10.**M_min_fit, 10.**M_1_fit, alpha_fit)
    #Nfit = Nall - Ndata[ibin_num]
    chi2 = sum(Fitsdata**2/ms_var**2) #+ Nfit**2/Nerr[ibin_num]**2
    print(M_min_fit, M_1_fit, alpha_fit)
    print('---chi2: %.3f'%chi2)
    return -0.5*chi2


#----------------------------#
def log_prior(theta):
    M_min_fit, M_1_fit, alpha_fit = theta
    if 10. < M_min_fit < 13. and 11. < M_1_fit < 15. and 0.5 < alpha_fit < 2.:
        return 0.
    return -np.inf


#----------------------------#
def log_probability(theta):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta)


#----------------------------#
print("---You're using the <Emcee_Fit> code---")
print("----------------------------------------")
 
#print("---Initial likelihood estimates:")
#Best_fits = np.loadtxt('./Data_Output/Fit_res/MiniChi2_HOD_' + bin_lum[ibin_num] + '.dat')
ini_x = np.array([10., 11., 0.5])
print(ini_x)

#----------------------------#
import emcee
num = 200
pos = np.zeros([num, 3])
pos[:, 0] = ini_x[0] + (3.*np.random.rand(num))
pos[:, 1] = ini_x[1] + (4.*np.random.rand(num))
pos[:, 2] = ini_x[2] + (1.5*np.random.rand(num))
pool = Pool()
nwalkers, ndim = pos.shape

if Paraller:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args = (), pool = pool)
else:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args = ())

sampler.run_mcmc(pos, 2000, progress = True)
#tau = sampler.get_autocorr_time()

labels = ["M_min", "M_1", "alpha_hod"]
flat_samples = sampler.get_chain(discard = 10, thin = 15, flat=True)
np.save(out_path + 'Fit_res/HOD_' + bin_lum[0] + '_' + bin_lum[ibin_num + 1] + '_samples_3param', flat_samples)


#----------------------------#
fig = corner.corner(flat_samples, labels=labels);
plt.savefig(out_path + 'Fit_res/HOD_' + bin_lum[0] + '_' + bin_lum[ibin_num + 1] + '_3param.png')

print(flat_samples)
from IPython.display import display, Math
#for i in range(ndim):
#    mcmc[i, :3] = np.percentile(flat_samples[i, :], [16, 50, 84])
#    mcmc[i, 3:] = np.diff(mcmc[i, :3])
txt_store = []
for i in range(ndim):
    mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
    q = np.diff(mcmc)
    txt = "\mathrm{{{3}}} = {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}"
    txt = txt.format(mcmc[1], q[0], q[1], labels[i])
    print(txt)
    display(Math(txt))
    txt_store.append(txt)

np.save(out_path + 'Fit_res/HOD_' + bin_lum[0] + '_' + bin_lum[ibin_num + 1] + '_3param', txt_store)


#----------------------------#
fig, axes = plt.subplots(ndim, figsize=(10, 7), sharex=True)
for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)
axes[-1].set_xlabel("step number");
plt.savefig(out_path + 'Fit_res/HOD_' + bin_lum[0] + '_' + bin_lum[ibin_num + 1] + '_3param_run_process.png')
plt.show()


print('---The end for <Emcee_Fit> code---')
print(' ')


