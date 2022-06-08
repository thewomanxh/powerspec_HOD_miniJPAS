import numpy as np
import sys
import getopt

## Here saves HOD parameters which need to fit.


#-----------Parameters initial value--------------#
M_min = 2.793e12#[1.e11, 1.e14]
## Unit: M_sun
alpha = 1.243#[0., 2.]
M_1 = 4.81e13#[1.e11.5, 1.e15.5]
## Unit: M_sun


#-------------JPAS and global factors-------------#
M_max_dm = 18.
## Unit: log10(M_sun)
M_min_dm = 10.
## Unit: log10(M_sun)

minzz = 0.01
## minimum redshift
maxzz = 2.
## maximum redshift
fsky = [1./(4.*np.pi)]
catntot = 7365.
## JPAS total catalog numbers
#bin_lum = ['-25.2', '-17.5', '-20.0', '-21.0', '-22.0', '-23.0']
bin_lum = ['-25.2', '-20.0', '-20.5', '-21.0', '-21.5', '-22.0', '-22.5', '-23.0']

#-------------------------------------------------#
verbose = False
## do you want print more information?
overwrite = False
Paraller = False
## do you want paraller computer when fitting?

redshift = np.linspace(minzz, maxzz, num = 30, endpoint = True)
wavenumber = np.linspace(1, 1000, num = 100, endpoint = True)
multipole = ([30, 47, 74, 117, 182, 284, 444, 692, 1078, 1680])
out_path = './Data_Output/'

"""
#--------------------ARGUMENTS-------------------#
opts, args = getopt.getopt(sys.argv[1:], \
            '-h-v-C:-M:-A:-O-P', \
            ['help', 'verbose', 'Mmin', 'M1', 'Alpha', 'overwrite', 'paraller'])

for opt_name, opt_value in opts:
    if opt_name in ('-h', '--help'):
        print(' ')
        print('---Welcome to use powespec_HOD_miniJPAS!')
        print(' -h --help           |you can read README at first.')
        print(' -v --verbose        |print more information. [default: False]')
        print(' -C --Mmin           |a parameter min Mass, format: *e**. [default:%.4e'%M_min, ']')
        print(' -M --M1             |a parameter Mass1, format: *e**. [default:%.4e'%M_1, ']')
        print(' -A --Alpha          |a parameter alpha. [default:', alpha, ']')
        print(' -O --Overwrite      |overwrite the output file.  [default: False]')
        print(' -P --paraller       |use Paraller. [default: False]')
        exit()

    if opt_name in ('-v', '--verbose'):
        verbose = True
    if opt_name in ('-C', '--Mmin'):
        M_min = float(opt_value)
    if opt_name in ('-M', '--M1'):
        M_1 = float(opt_value)
        alpha = float(opt_value)
    if opt_name in ('-O', '--overwrite'):
        overwrite = True
    if opt_name in ('-P', '--paraller'):
        Paraller = True        
"""

#---------------------PRINT----------------------#
print(' ')
#print('---The initial parameters are in JPAS:')
#print('   |M_min: %.3e'%M_min)
#print('   |M_1: %.3e'%M_1)
#print('   |alpha: %.3f'%alpha)
#print(' ')


