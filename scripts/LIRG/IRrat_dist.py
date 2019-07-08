"""
Emily Biermann
Ratio of bright IR objects to blue cloud objects as function of distance from
cluster center.
7/28/19
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import stats
from scipy.stats import norm

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

# Save Figures?
save = False
if save:
    tag = 'LIRG' # tag for figure names

# Show Figures?
show = True

# Data Directories

# IR Bright
pwd1  = '/home/ebiermann/cat/LIRG/'
fname1 = 'ms0451_crawSpec_MACS0454_match.fits'

# Blue Cloud (BC)
pwd2 = '/home/ebiermann/cat/MACS0454-3_typed/'
fname2 = 'WtGMACS0454cut-SpecCraw_Selection.fits'

# All
pwd3 = '/home/ebiermann/cat/MACS0454-3_typed/'
fname3 = 'galaxySpecCraw_ClusterZ.fits'

cat1 = fits.open(pwd1 + fname1)
data1 = cat1[1].data
cat1.close()

cat2 = fits.open(pwd2 + fname2)
data2 = cat2[1].data
cat2.close()

cat3 = fits.open(pwd3 + fname3)
data3 = cat3[1].data
cat3.close()

d_IR = data1['RADIUS']
d_BC = data2['RADIUS']
d_all = data3['RADIUS']

# Find Ratios at step intervals

rat = []
dist_mid = []

dmax = np.amax(d_all)
dmin = np.amin(d_all)
step = 40

lower = dmin
while lower < dmax:
    upper = lower+step
    objs = np.where(np.logical_and(d_IR>=lower,d_IR<upper))
    num_IR = len(objs[0])
    objs = np.where(np.logical_and(d_BC>=lower,d_BC<upper))
    num_BC = len(objs[0])
    if num_BC != 0:
        rat.append(float(num_IR)/float(num_BC))
        dist_mid.append((lower+upper)/2.0)
    lower = lower + step

# Plot ratio vs. distance

plt.figure()
plt.title('Galaxy Cluster MACS0454')
plt.xlabel('Distance From Cluster Center, interval={} (arcsec)'.format(step))
plt.ylabel('Ratio of IR-Bright to Blue Cloud')
plt.errorbar(dist_mid,rat,fmt='.')
if save:
    plt.savefig('../figures/LIRG/ratVsDist_{}.png'.format(tag),format='png',dpi=1000,bbox_inches='tight')
    

# Show Figure
if show:
    plt.show()
plt.clf

