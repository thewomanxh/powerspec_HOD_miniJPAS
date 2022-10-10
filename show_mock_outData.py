import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
import h5py
import cartopy.crs as ccrs

from plotParam import *


## This code is to show the output galaxy catalogue mocks.

fpath = './Data_Output/'
fname = 'galaxy_catalogue_pinocchio.hdf5'
show_comp = True
show_hist = True


#----------------------------------#
def convert_ra2long(ravals):
    ## to convert RA values defined in [0, 360] to LONG values defined in [-180, 180]
    longvals = np.atleast_1d(ravals.copy())
    longvals[longvals > 180.] = longvals[longvals > 180] - 360.
    return longvals

outFile = './hodpy_pinocchio/output/' + fname
outData = h5py.File(outFile)
inFile = './hodpy_pinocchio/input/halo_catalogue_pinocchio.hdf5'
inData = h5py.File(inFile)


#-------------Plot-----------------#
if show_comp:
    f = plt.figure(1, (38, 36))
    ax = f.add_subplot(111, projection = ccrs.AzimuthalEquidistant(central_latitude = 90.))
    ax.gridlines(draw_labels = True, y_inline = True, xformatter = ticker.StrMethodFormatter("{x:.3g}"), yformatter = ticker.StrMethodFormatter("{x:.3g}"))
    ax.scatter(convert_ra2long(np.array(inData['Data']['ra'][::2000])), inData['Data']['dec'][::2000], s=2, transform=ccrs.PlateCarree(), label = 'Input Halo Catalogue', marker = '.', alpha = 0.5, rasterized = True)
    ax.scatter(convert_ra2long(np.array(outData['ra'][::2000])), outData['dec'][::2000], s=1, transform=ccrs.PlateCarree(), label = 'Output Galaxy Catalogue', c = 'red', marker = 'o', alpha = 0.5, rasterized=True) 
    plt.savefig(fpath + 'galaxy_catalogue_pinocchio_compare.png')
    ax.legend(loc = 'lower left', bbox_to_anchor = (0, 1), frameon = False)
    plt.show()


if show_hist:
    f, ax = plt.subplots(figsize = (20, 18))
    ax.hist(outData['zcos'][::3000], bins = 100, density = True, label = 'Output Galxy Catalogue')
    plt.xlim([0.2, 0.6])
    plt.xticks([0.2, 0.3, 0.4, 0.5, 0.6], ('0.2', '0.3', '0.4', '0.5', '0.6'), fontsize = 12)
    ax.set_xlabel(r'$z$', fontsize = 15)
    ax.set_ylabel('Numbers', fontsize = 15) 
    ax.legend(loc = 'upper right', frameon = False)
    f.savefig(fpath + 'galaxy_catalogue_pinocchio_distribution.png')
    plt.show()


