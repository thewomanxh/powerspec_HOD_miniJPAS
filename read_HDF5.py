import sys
import os
import getopt
import numpy as np
import h5py
##Here is the code to read and show data form HDF5 files.


#-----------------ARGUMENTS-----------------#
opts, args = getopt.getopt(sys.argv[1:], \
             '-h-p:-f:', \
             ['help', 'path', 'file-name'])

for opt_name, opt_value in opts:
    if opt_name in ('-h', '--help'):
        print(' ')
        print(' -h --help')
        print(' -p --path           |input file path.')
        print(' -f --file-name      |input file name.')

        exit()
    if opt_name in ('-p', '--path'):
        fpath = str(opt_value)
    if opt_name in ('-f', '--file-name'):
        fname = str(opt_value)


#-----------------OPERATION-----------------#
#fname = './input/halo_catalogue_small.hdf5'
fdata = h5py.File(fpath + fname, 'r')
opath = './Data_Output/'

print('-------------------------------------')
print(' The keys in file ', fname, 'are (is):')
print(' ',list(fdata.keys()))
print('-------------------------------------')

big_idex = 1
small_idex = 1

if '/input/' in fpath:
    for group in fdata.keys():
        print('{', big_idex, '}', group)
        print(' ')
        
        for dset in fdata[group].keys():
            print('(', small_idex,')', dset)
            ds_data = fdata[group][dset]
            ## returns HDF5 dataset object
            print(ds_data)
            #print(ds_data.shape, ds_data.dtype)

            arr = fdata[group][dset][:]
            ## adding [:] returns a numpy array
            if dset == 'dec':
                np.savetxt(opath + 'halo_dec.txt', arr, '%3f')
            if dset == 'ra':
                np.savetxt(opath + 'halo_ra.txt', arr, '%3f')
 
            print(arr.shape, arr.dtype)
            print(arr)
            print(' ')
            small_idex = small_idex + 1
    big_idex = big_idex + 1

else:
    with open(opath + fname + ".txt", "ab") as ff:
        for group in fdata.keys():
            print('{', big_idex, '}', group)
            print(' ')
            print('(', small_idex,')', group)
            ds_data = fdata[group]
            print(ds_data)

            arr = fdata[group][:]
            print(arr.shape, arr.dtype)
            print(arr)
            print(' ')
            if small_idex == 1:
                res = np.zeros(arr.shape[0])
            res = np.column_stack((res, arr))
            small_idex = small_idex + 1
            
            if group == 'is_cen':
                cen = np.count_nonzero(arr)
        print('---total central galaxy is', cen)

        res = np.delete(res, 0, axis = 1)
        head = ",".join(fdata.keys())
        np.savetxt(ff, res, fmt = '%.6e', header=head) 
        big_idex = big_idex + 1



