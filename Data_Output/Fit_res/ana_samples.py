import os
import numpy as np
import time
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from multiprocessing import Pool
import corner
import sys
os.environ["OMP_NUM_THREADS"] = "1" 


#-----Global parameters------#
global bin_lum, ibin_num
bin_lum = ['-25.2', '-20.0', '-20.5', '-21.0', '-21.5', '-22.0', '-22.5', '-23.0']
#bin_lum = ['-25.2', '-17.5', '-20.0', '-21.0', '-22.0', '-23.0']
ibin_num = int(sys.argv[1])

samples=np.load('./HOD_' + bin_lum[0] + '_' + bin_lum[ibin_num + 1] + '_samples_3param.npy')

nwalkers = 200
nrun = 2000
labels = ["M_min", "M_1", "alpha_hod"]
flat_samples = (samples.reshape((200*2000,3)))
#flat_samples = samples.flatten(1)#.get_chain(discard = 10, thin = 15, flat=True)
ndim = 3

#----------------------------#
fig = corner.corner(flat_samples, labels=labels);
plt.savefig('./HOD_' + bin_lum[0] + '_' + bin_lum[ibin_num + 1] + '_3param.png')

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

np.save('./HOD_' + bin_lum[0] + '_' + bin_lum[ibin_num + 1] + '_3param', txt_store)


#----------------------------#
fig, axes = plt.subplots(ndim, figsize=(10, 7), sharex=True)
for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)
axes[-1].set_xlabel("step number");
plt.savefig('./HOD_' + bin_lum[0] + '_' + bin_lum[ibin_num + 1] + '_3param_run_process.png')
#plt.show()


print('---The end for <Emcee_Fit> code---')
print(' ')


