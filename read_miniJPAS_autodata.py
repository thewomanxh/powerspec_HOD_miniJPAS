import numpy as np
import pandas as pd
## This code is to read miniJPAS auto correlation calculated by Luis.

fpath = './Data_Input/JPAS/'
fname = False#'miniJPAS_BB_LSS.auto.txt'
fcsv = 'minijpas_pdr201912_cat_wAbsMag_v1.csv'
fout = './Data_Output/miniJPAS_selection/'

band_Switch = [8, 10, 12, 14, 16, 18]
band_name = ['uJAVA', 'uJPAS', 'g', 'r', 'i', 'J10066']
## Order is: uJAVA, uJPAS, g, r, i, J10069


#-------------DATA---------------#
#with open(fpath + fname) as f:
#    for line in f.readlines():
#        print(line[14])
if fname:
    fdata = np.loadtxt(fpath + fname)
    Del_ind = []
    for i in range(len(fdata)):
        if (fdata[i, 6] > 0.1 and fdata[i, 7] > 0.1):
            Del_ind.append(i)
    fdata_no_star = np.delete(fdata, Del_ind, axis = 0)
    
    Mr = [min(fdata_no_star[:, 16]), -23., -22., -21., -20., max(fdata_no_star[:, 16])]
    print(Mr)
    fdata_i_index = np.zeros([len(fdata), len(Mr)-1])
    
    for i in range(len(Mr)-1):
        for j in range(len(fdata_no_star)):
            if Mr[0] > fdata_no_star[j, 16] or fdata_no_star[j, 16] > Mr[i+1]:
                fdata_i_index[j, i] = j
        f_i = np.delete(fdata_no_star, np.nonzero(fdata_i_index[:, i]), axis = 0)
        print(i, np.nonzero(fdata_i_index[:, i]),f_i, len(f_i))
        np.savetxt(fout + 'miniJPAS_BB_LSS_iband_selection_%.1f'%Mr[i] + '.txt', f_i, \
                   fmt = '%d  %d  %d  %d  %8f   %8f  %3f  %3f  %4f  %5f  %4f  %5f  %4f  %5f  %4f  %5f  %4f  %5f  %4f  %5f  %1f')
else:
    fdata = np.array(pd.read_csv(fpath + fcsv))
    zz = fdata[:, 33]
    ## remain redshift from 0.2-0.6
    Del_zz_ind = []; Del_star_ind = []
    for i in range(len(fdata)):
        if (fdata[i, 7] > 0.5):
            Del_star_ind.append(i)
        if (fdata[i, 33] > 0.6 or fdata[i, 33] < 0.2):
            Del_zz_ind.append(i)
    Del_ind = list(set(Del_zz_ind).union(set(Del_star_ind)))
    fdata_zz_nostar = np.delete(fdata, Del_ind, axis = 0)

    Mr = [min(fdata_zz_nostar[:, 16]), -23., -22., -21., -20., max(fdata_zz_nostar[:, 16])]
    fdata_i_index = np.zeros([len(fdata), len(Mr)])

    for i in range(len(Mr)-1):
        for j in range(len(fdata_zz_nostar)):
            if Mr[0] > fdata_zz_nostar[j, 16] or fdata_zz_nostar[j, 16] > Mr[i+1]:
                fdata_i_index[j, i] = j
        f_i = np.delete(fdata_zz_nostar, np.nonzero(fdata_i_index[:, i]), axis = 0)
        f_i_j = np.delete(f_i, [1, 2, 3, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50], axis = 1)
        np.savetxt(fout + 'miniJPAS_zz_nostar_iband_selection_%.1f'%Mr[0] + '_%.1f'%Mr[i+1] + '.txt', f_i_j, fmt = '%d    %8f    %8f    %5f    %8f     %8f     %8f')




