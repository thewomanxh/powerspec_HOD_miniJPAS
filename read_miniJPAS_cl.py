import numpy as np

## This code is to read miniJPAS Cl calculated by Carlos.

fpath = './Data_Input/Measurement_Cl_from_carlos/'
fname_1 = '/msk_miniJPAS_dth_1arcmin/'
sname_1 = ['map_msk', 'phiR0', 'thetaR0', 'Rmatrix']
fname_2 = '/pofk_miniJPAS_v2_finepix/'
sname_2 = ['map_msk', 'k1d', 'pk_map', 'pk_Poisson', 'dmap_out']
name_opt_1 = [1, 1, 1, 1]
name_opt_2 = [0, 0, 0, 0, 1]
lum_out = ['-25.2', '-23.0', '-22.0', '-21.0', '-20.0', '-17.5']
fout = './Data_Output/miniJPAS_cl/'
cl_opt = 1
error_opt = 1

#-----------READ DATA------------#
for i in range(len(name_opt_1)):
    if name_opt_1[i] == 1:
        fdata = np.load(fpath + fname_1 + sname_1[i] + '.npy', encoding = 'latin1')

for j in range(len(name_opt_2)):
    if name_opt_2[j] == 1:
        fdata = np.load(fpath + fname_2 + sname_2[j] + '.npy', encoding = 'latin1')

k1d = np.load(fpath + fname_2 + sname_2[1] + '.npy', encoding = 'latin1')
pk_map = np.load(fpath + fname_2 + sname_2[2] + '.npy', encoding = 'latin1')
pk_poisson = np.load(fpath + fname_2 + sname_2[3] + '.npy', encoding = 'latin1')


#----------------Cl--------------#
if cl_opt == 1:
    effictive_cl = 0*pk_map
    for k in range(5):
        effictive_cl[k, :] = pk_map[k, :] - np.mean(pk_poisson[k, :, :], axis = 0)
        np.savetxt(fout + 'miniJPAS_' + lum_out[0] + '_' + lum_out[k+1] + '_cl.dat', np.column_stack((k1d, effictive_cl[k, :])))
    

#-------------Error--------------#
if error_opt == 1:
    for k in range(5):
        relative_error = np.std(pk_poisson[k, :, :], axis = 0) / np.mean(pk_poisson[k, :, :])
        relative_error_obs = pk_map[k, :] *relative_error
        np.savetxt(fout + 'miniJPAS_' + lum_out[0] + '_' + lum_out[k+1] + '_error.dat', relative_error_obs)


