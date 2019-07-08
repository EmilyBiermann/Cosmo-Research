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

pwd  = '/home/ebiermann/cat/MACS0454-3_typed/'
fname_match = 'galaxySpecCraw_ClusterZ.fits'
fname_red = 'galaxySpecCraw_ClusterZ_redSeq.fits'
fname_green = 'galaxySpecCraw_ClusterZ_greenVal.fits'
fname_blue = 'galaxySpecCraw_ClusterZ_blueCloud.fits'

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

# All objects

matchCat = fits.open(pwd + fname_match)
data = matchCat[1].data
matchCat.close()

c = 2.988E8 # m/s

B = data['MAG_APER1-SUBARU-conv-1-W-J-B']
V = data['MAG_APER1-SUBARU-conv-1-W-J-V']
I = data['MAG_ISO-SUBARU-conv-1-W-S-I+']

color = B - V

d_center = data['RADIUS']
velocity = c*data['Z']

plt.figure()
plt.xlabel('I')
plt.ylabel('B-V')
plt.errorbar(I,color,fmt='.')

#-------------------------------------------------------------------------------

# By color

redCat = fits.open(pwd + fname_red)
data = redCat[1].data
redCat.close()
B_r = data['MAG_APER1-SUBARU-conv-1-W-J-B']
V_r = data['MAG_APER1-SUBARU-conv-1-W-J-V']
I_r = data['MAG_ISO-SUBARU-conv-1-W-S-I+']
color_r = B_r - V_r
d_r = data['RADIUS']
vel_r = c*data['Z']

greenCat = fits.open(pwd + fname_green)
data = greenCat[1].data
greenCat.close()
B_g = data['MAG_APER1-SUBARU-conv-1-W-J-B']
V_g = data['MAG_APER1-SUBARU-conv-1-W-J-V']
I_g = data['MAG_ISO-SUBARU-conv-1-W-S-I+']
color_g = B_g - V_g
d_g = data['RADIUS']
vel_g = c*data['Z']

blueCat = fits.open(pwd + fname_blue)
data = blueCat[1].data
blueCat.close()
B_b = data['MAG_APER1-SUBARU-conv-1-W-J-B']
V_b = data['MAG_APER1-SUBARU-conv-1-W-J-V']
I_b = data['MAG_ISO-SUBARU-conv-1-W-S-I+']
color_b = B_b - V_b
d_b = data['RADIUS']
vel_b = c*data['Z']

plt.figure()
plt.xlabel('I')
plt.ylabel('B - V')
plt.errorbar(I_r,color_r,fmt='.',color='red',label='Red Sequence')
plt.errorbar(I_g,color_g,fmt='.',color='green',label='Green Valley')
plt.errorbar(I_b,color_b,fmt='.',color='blue',label='Blue Cloud')
plt.legend()

plt.figure()
plt.xlabel('Distance From Cluster Center (arcsec)')
plt.ylabel('Velocity')
plt.errorbar(d_r,vel_r,fmt='.',color='red',label='Red Sequence')
plt.errorbar(d_g,vel_g,fmt='.',color='green',label='Green Valley')
plt.errorbar(d_b,vel_b,fmt='.',color='blue',label='Blue Cloud')
plt.legend()

if show:
    plt.show()
plt.clf()
