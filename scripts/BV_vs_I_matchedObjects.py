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
save = True
if save:
    tag = 'MACS0454_13match' # tag for figure names

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
# Topcat match MUST have WtG as first catalog!!

pwd  = "/home/ebiermann/catalogs/"
#fname_match = 'WtG_SpecCrawNoStar_1arc_3arc_match.fits'
#cat1=True
#cat2=False

fname_match = 'MACS0454-3_typed/WtGMACS0454-03typed_SpecCrawNoStar_13match.fits'
cat1=False
cat2=True

#fname_match = 'WtG_SpecCrawNoStar_3arcsecMatch.fits' #Craw: ZTYPE==1, STAR==0
#fname_match = 'WtGPhot_CSpec_match.fits' # 1arcsec, stars included
#fname_match = 'WtG_SpecCrawNoStarNonzeroZ_3arcsecMatch.fits'
#fname_match = 'WtG_SpecCrawNoStarNonzeroZ_1arcsecMatch.fits'
#fname_match = 'WtG_SpecCraw_magStarZcuts_1arcsecMatch.fits'
#fname_match = 'WtG_SpecCraw_magStarZcuts_3arcsecMatch.fits'

matchCat = fits.open(pwd + fname_match)
data = matchCat[1].data

# Allocate arrays for magnitude, color, redshift
WtG_BV = []
WtG_I = []
specZ = []
specZ_cut = []

# Get data
if cat1:
    for i in range(0,len(data)):
        z = data[i][774] # Z (775)
        specZ.append(z)
        if z>=0.524 and z<=0.552:
            WtG_B = data[i][528]  # MAG_APER1_SUBARU-10_2-1-W-J-B (529)
            WtG_V = data[i][556]  # MAG_APER1_SUBARU-10_2-1-W-J-V (557)
            WtG_BV.append(WtG_B - WtG_V)
            WtG_I.append(data[i][616])  # MAG_ISO-SUBARU-10_2-1-W-S-I+ (617)
            specZ_cut.append(z)
if cat2:
    for i in range(0,len(data)):
        z = data[i][829] # Z (830)
        specZ.append(z)
        #if z>=0.524 and z<=0.552:
        if cat2:
            WtG_B = data[i][603]  # MAG_APER1_SUBARU-conv-1-W-J-B (604)
            WtG_V = data[i][623]  # MAG_APER1_SUBARU-conv-1-W-J-V (624)
            WtG_BV.append(WtG_B - WtG_V)
            WtG_I.append(data[i][675])  # MAG_ISO-SUBARU-conv-W-S-I+ (676)
            specZ_cut.append(z)

# Histogram
nbins = 50

# SpecZ Histogram to determine cluster redshift
mean = np.average(specZ_cut)
stdev = np.std(specZ_cut)
print('Unclipped')
print('mean =', mean)
print('stdev = ', stdev)
print(' ')
xarray,yarray,mode,cmin,cmax =\
     gaussFit(specZ_cut,nbins,mean,stdev)
# Plot Figure
plt.figure()
plt.title(r'SpecZ')
plt.xlabel(r'Spectroscopic Redshift')
plt.ylabel(r'Number of Galaxies')
plt.hist(specZ_cut,range=[cmin,cmax], bins=nbins);
plt.plot(xarray,yarray,color="black",linewidth=1.0,label='Gaussian Fit')
plt.axvline(x=mean,linewidth=1.0,color="yellow",label='mean')
plt.axvline(x=mean-stdev,linewidth=1.0,color='red',label='stdev')
plt.axvline(x=mean+stdev,linewidth=1.0,color='red')
plt.legend()
if save:
    plt.savefig('../figures/SpecZhist_{}.png'.format(tag),format='png',dpi=1000,bbox_inches='tight')
'''
# Color Magnitude Diagram
plt.figure()
plt.title(r'Color Magnitude Comparison')
plt.xlabel(r'I')
plt.ylabel(r'B - V')
plt.errorbar(WtG_I, WtG_BV, fmt='.')
if save:
    plt.savefig('../figures/BVvsI_{}.png'.format(tag),format='png',dpi=1000,bbox_inches='tight')
'''
if show:
    plt.show()
plt.clf
