import numpy as np
import sys
import getopt

## Here saves HOD parameters which need to fit.


#-----------Parameters initial value--------------#
M_min = 2.793e12#[1.e11, 1.e14]
## Unit: M_sun
alpha = 1.243#[0., 2.]
M_1 = 4.81e13#[1.e11.5, 1.e15.5]
## Unit: M_sun, if it's the same with M_0, the model is 4 parameters model


#-------------JPAS and global factors-------------#
M_max = 16.
## Unit: log10(M_sun)
M_min = 10.
## Unit: log10(M_sun)

minzz = 0.01
maxzz = 2.
fsky = [1./(4.*np.pi)]
catntot = 7365.
## JPAS total catalog numbers.


#-------------------------------------------------#
Fit = False
verbose = False
overwrite = False
Paraller = True

redshift = np.linspace(0., 1.5, num = 30, endpoint = True)
wavenumber = np.linspace(1, 1000, num = 100, endpoint = True)
multipole = ([30, 47, 74, 117, 182, 284, 444, 692, 1078, 1680])
out_path = './Data_Output/'


#--------------------ARGUMENTS-------------------#
opts, args = getopt.getopt(sys.argv[1:], \
            '-h-v-C:-M:-A:-O-F-P', \
            ['help', 'verbose', 'Mmin', 'M1', 'Alpha', 'overwrite', 'Fit', 'paraller'])

for opt_name, opt_value in opts:
    if opt_name in ('-h', '--help'):
        print(' ')
        print(' -h --help')
        print(' -v --verbose        |[default: True]')
        print(' -C --Mmin           |a parameter min Mass, format: *e**. [default:', M_min, ']')
        print(' -M --M1             |a parameter Mass1, format: *e**. [default:', M_1, ']')
        print(' -A --Alpha          |a parameter alpha. [default:', alphs, ']')
        print(' -O --Overwrite      |Overwrite the output file.  [default: False]')
        print(' -F --Fit            |use Emcee MCMC. [default: False]')
        print(' -P --paraller       |use Paraller. [default: False]')
        exit()

    if opt_name in ('-v', '--verbose'):
        verbose = True
    if opt_name in ('-C', '--Mmin'):
        M_min = float(opt_value)
    if opt_name in ('-M', '--M1'):
        M_1 = float(opt_value)
        alphs = float(opt_value)
    if opt_name in ('-O', '--overwrite'):
        overwrite = True
    if opt_name in ('-F', '--Fit'):
        Fit = True
    if opt_name in ('-P', '--paraller'):
        Paraller = True        


#---------------------PRINT----------------------#
print(' ')
#print('---The initial parameters are in JPAS:')
#print('   |M_min: %.3e'%M_min)
#print('   |M_1: %.3e'%M_1)
#print('   |alphs: %.3f'%alphs)
#print(' ')

