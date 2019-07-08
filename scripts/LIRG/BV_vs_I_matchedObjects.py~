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

# Save Figures?
save = False
if save:
    tag = 'LIRG' # tag for figure names

# Show Figures?
show = True

# Definitions
def gaussFit(data,nbins,mu,sig):
    cmin=np.amin(data)
    cmax=np.amax(data)
    nbins=100
    mode=stats.mode(data)[0][0]
    normalization=(cmax-cmin)/nbins*len(data)
    xarray=np.linspace(cmin,cmax,nbins*10)
    yarray=normalization*norm.pdf(xarray,loc=mu, scale=sig)
    return(xarray,yarray,mode,cmin,cmax)

# Data Directories

pwd1  = "/home/ebiermann/cat/LIRG/"

fname_match = 'ms0451_crawSpec_MACS0454_match.fits'
Bname = 'MAG_APER1-SUBARU-conv-1-W-J-B'
Vname = 'MAG_APER1-SUBARU-conv-1-W-J-V'
Iname = 'MAG_ISO-SUBARU-conv-1-W-S-I+'

pwd2 = '/home/ebiermann/cat/MACS0454-3_typed/'
fname2 = 'galaxySpecCraw_ClusterZ.fits'

matchCat = fits.open(pwd1 + fname_match)
data = matchCat[1].data
matchCat.close()

cat2 = fits.open(pwd2 + fname2)
data2 = cat2[1].data
cat2.close()

# Allocate arrays for magnitude, color, redshift
BV_1 = data[Bname] - data[Vname]
I_1 = data[Iname]
specZ_1 = data['Z']

BV_2 = data2[Bname] - data2[Vname]
I_2 = data2[Iname]
specZ_2 = data2['Z']


# Histogram
nbins = 50
'''
# SpecZ Histogram to determine cluster redshift
mean = np.average(specZ)
stdev = np.std(specZ)
print('Unclipped')
print('mean =', mean)
print('stdev = ', stdev)
print(' ')
xarray,yarray,mode,cmin,cmax =\
     gaussFit(specZ,nbins,mean,stdev)
# Plot Figure
plt.figure()
plt.title(r'SpecZ')
plt.xlabel(r'Spectroscopic Redshift')
plt.ylabel(r'Number of Galaxies')
plt.hist(specZ,range=[cmin,cmax], bins=nbins);
plt.plot(xarray,yarray,color="black",linewidth=1.0,label='Gaussian Fit')
plt.axvline(x=mean,linewidth=1.0,color="yellow",label='mean')
plt.axvline(x=mean-stdev,linewidth=1.0,color='red',label='stdev')
plt.axvline(x=mean+stdev,linewidth=1.0,color='red')
plt.legend()
if save:
    if cat1:
        plt.savefig('../figures/SpecZhist_{}.png'.format(tag),format='png',dpi=1000,bbox_inches='tight')
    if cat2:
        plt.savefig('../figures/MACS0454-03_typed/SpecZhist_{}.png'.format(tag),format='png',dpi=1000,bbox_inches='tight')
'''
# Color Magnitude Diagram
plt.figure()
plt.title(r'Color Magnitude Comparison')
plt.xlabel(r'I')
plt.ylabel(r'B - V')
plt.xlim(17.0,24.6)
plt.ylim(0.,1.8)
plt.errorbar(I_2, BV_2, fmt='.')
plt.errorbar(I_1, BV_1, fmt='.',label='FIR Objects')
plt.legend()
if save:
    plt.savefig('../figures/BVvsI_{}.png'.format(tag),format='png',dpi=1000,bbox_inches='tight')

if show:
    plt.show()
plt.clf
