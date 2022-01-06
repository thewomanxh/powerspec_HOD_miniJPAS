import matplotlib as mp
import matplotlib.pyplot as plt 
from matplotlib.ticker import FormatStrFormatter
## This file is to the parameter for xiu's plots.

mp.rcParams['lines.linewidth'] = 2.
mp.rcParams['ytick.right']=True
mp.rcParams['xtick.top']=True
mp.rcParams['ytick.direction']='in'
mp.rcParams['xtick.direction']='in'
mp.rcParams['ytick.major.width']=1
mp.rcParams['xtick.major.width']=1
mp.rcParams['ytick.minor.width']=1
mp.rcParams['xtick.minor.width']=1
mp.rcParams['ytick.minor.size']=4.
mp.rcParams['xtick.minor.size']=4.
mp.rcParams['xtick.major.pad']=10.


ax = plt.gca()
ax.xaxis.get_minor_ticks()
ax.spines['bottom'].set_linewidth(1)
ax.spines['left'].set_linewidth(1)
ax.spines['right'].set_linewidth(1)
ax.spines['top'].set_linewidth(1)
ax.xaxis.set_major_formatter(FormatStrFormatter('%i'))
#plt.legend(loc = 'best', fontsize = 30, frameon = False)
#print mp.rcParams

