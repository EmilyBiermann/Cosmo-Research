"""
Emily Biermann
Compare mag from WtG, Crawford Match
3/14/19
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

#-------------------------------------------------------------------------------

cat1 = True
cat2 = False

# Data Directories
# Topcat match MUST have WtG as first catalog!!

pwd  = "/home/ebiermann/cat/mag_matchCat/"
if cat1:
    fname_match = 'match4_1as1mag_3as1mag.fits'
    
if cat2:
    fname_match = 'match3_3as1mag_rem.fits'

# Save Figures?
save = False
'''
if save:
    if cat1:
        pwdfig = '/home/ebiermann/Cosmo-Research/figures/mag_matchCat/all/galaxy/'
        tag = 'all' # tag for figure names
    if cat2:
        pwdfig = '/home/ebiermann/Cosmo-Research/figures/mag_matchCat/match3_rem/mag/'
        tag = 'match3rem'
'''

# Show Figures?
show = True

#-------------------------------------------------------------------------------

matchCat = fits.open(pwd + fname_match)
data = matchCat[1].data

d_center = []
velocity = []

c = 2.988E8 # m/s

for i in range(0,len(data)):
    z = data[i][776] # Z, 777
    if z>=0.524 and z<=0.552:
        d_center.append(data[i][679]) # RADIUS, 680
        velocity.append(c*z)

plt.figure()
plt.xlabel('Distance From Cluster Center (arcsec)')
plt.ylabel('Velocity')
plt.errorbar(d_center,velocity,fmt='o')
if show:
    plt.show()
plt.clf()
