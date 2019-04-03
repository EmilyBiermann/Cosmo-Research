"""
Emily Biermann
Photo vs Spec z
3/4/19
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
    tag = 'final' # tag for figure name

# Show Figures?
show = True

# Defintions

def sigma3clip(data,mean,stdev):
    # Clips data over 3sigma from mean, 
    # returns aray without these points, number of data points cut
    # returns new mean, stdev
    dataClip = []
    numCut = 0
    for d in data:
        if np.abs(d-mean)<3*stdev:
            dataClip.append(d)
        else:
            numCut = numCut + 1
    return(dataClip,numCut,np.mean(dataClip),np.std(dataClip))

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

pwd  = "../cat/"
fname_match = 'WtG_SpecCrawNoStar_1arc_3arc_match.fits'

#fname_match = 'WtG_SpecCrawNoStar_3arcsecMatch.fits' #Craw: ZTYPE==1, STAR==0
#fname_match = 'WtGPhot_CSpec_match.fits' # 1arcsec, stars included
#fname_match = 'WtG_SpecCraw_Zcut_3arcsecMatch.fits' # Craw: 0.5<=z<=0.9
#fname_match = 'WtG_SpecCrawNoStarNonzeroZ_3arcsecMatch.fits'
#fname_match = 'WtG_SpecCrawNoStarNonzeroZ_1arcsecMatch.fits'
#fname_match = 'WtG_SpecCraw_magStarZcuts_1arcsecMatch.fits'
#fname_match = 'WtG_SpecCraw_magStarZcuts_3arcsecMatch.fits'


matchCat = fits.open(pwd + fname_match)
data = matchCat[1].data

bpzCat = fits.open(pwd + 'MACS0454-03.W-C-RC.bpz.tab')
bpzData = bpzCat[1].data

specZ = []
specZ_err = []

photoZ = []
photoZ_errAbove = []
photoZ_errBelow = []

photoZ_chi2 = []

for i in range(0,len(data)):
    specZ.append(data[i][774])
    specZ_err.append(data[i][775])
    
    WtG_SeqNr = data[i][2]             # Get data number
    pZ = bpzData[WtG_SeqNr-1][1]       # photoZ (BPZ_Z_B)
    photoZ.append(pZ)                  # Add photoZ to array
    pZmin = bpzData[WtG_SeqNr-1][2]    # Get min photoZ value
    photoZ_errBelow.append(pZ - pZmin) # Calculate lower errbar
    pZmax = bpzData[WtG_SeqNr-1][3]    # Get max photoZ value
    photoZ_errAbove.append(pZmax - pZ) # Calculate upper errbar
    photoZ_chi2.append(bpzData[WtG_SeqNr-1][8])

photoZ_err = np.array([photoZ_errBelow,photoZ_errAbove])

# 1-1 line
x=np.linspace(0,np.amax(specZ))
y=x

plt.figure()
plt.title(r'Redshift Comparison')
plt.xlabel(r'SpecZ from Crawford')
plt.ylabel(r'PhotoZ from WtG')
#plt.errorbar(specZ, photoZ, xerr=specZ_err, yerr=photoZ_err, fmt='x')
plt.errorbar(specZ, photoZ, xerr=specZ_err, fmt='.')
plt.plot(x,y,linestyle='--')
if save:
    plt.savefig('SpecPhoto_line_{}.png'.format(tag),format='png',dpi=1000,bbox_inches='tight')

# Histogram
nbins = 100

# Unclipped Data
points = []
i=0
for pZ in photoZ:
    points.append(pZ - specZ[i])
    i=i+1
mean = np.average(points)
stdev = np.std(points)
xarray,yarray,mode,cmin,cmax =\
     gaussFit(points,nbins,mean,stdev)
print('Unclipped')
print('mean =', mean)
print('stdev = ', stdev)
print(' ')
# Plot Figure
plt.figure()
plt.title(r'Unclipped Data')
plt.xlabel(r'Photometric Redshift - Spectroscopic Redshift')
plt.ylabel(r'Number of Galaxies')
plt.hist(points,range=[cmin,cmax], bins=nbins);
plt.plot(xarray,yarray,color="red",linewidth=1.0,label='Gaussian Fit')
plt.axvline(x=mean,linewidth=1.0,color="yellow",label='mean')
plt.legend()
text_stats = r'''\noindent $\mu = {:.3f}$ \\'''.format(mean) + \
             r'''$\sigma = {:.3f}$ \\'''.format(stdev)
#plt.text(1.3,325,text_stats,{'fontsize':20})
if save:
    plt.savefig('SpecPhoto_noclip_{}.png'.format(tag),format='png',dpi=1000,bbox_inches='tight')

### Sigma Clipping 1 ###
points_clip1,pointsCut1,mean_clip1,stdev_clip1 = sigma3clip(points,mean,stdev)
totclip = pointsCut1
print('3sigma clipped, iter1:',pointsCut1)
print('mean =', mean_clip1)
print('stdev = ', stdev_clip1)
print(' ')
# Fit Gaussian
xarray,yarray,mode,cmin,cmax =\
     gaussFit(points_clip1,nbins,mean_clip1,stdev_clip1)
# Plot Figure
plt.figure()
plt.title(r'1 Iteration of 3 Sigma Clip')
plt.xlabel(r'Photometric Redshift - Spectroscopic Redshift')
plt.ylabel(r'Number of Galaxies')
plt.hist(points_clip1,range=[cmin,cmax], bins=nbins);
plt.plot(xarray,yarray,color="red",linewidth=1.0,label='Gaussian Fit')
plt.axvline(x=mean_clip1,linewidth=1.0,color="yellow",label='mean')
plt.legend()
text_stats = r'''\noindent $\mu = {:.3f}$ \\'''.format(mean_clip1) + \
             r'''$\sigma = {:.3f}$ \\'''.format(stdev_clip1)
text_clip = r'''Total Clippied = {}'''.format(totclip)
#plt.text(0.4,175,text_stats,{'fontsize':20})
#plt.text(-0.6,225,text_clip)
if save:
    plt.savefig('SpecPhoto_clip1_{}.png'.format(tag),format='png',dpi=1000,bbox_inches='tight')

### Sigma Clipping 2 ###
points_clip2,pointsCut2,mean_clip2,stdev_clip2 = \
    sigma3clip(points_clip1,mean_clip1,stdev_clip1)
totclip = totclip + pointsCut2
print('3sigma clipped, iter2:',pointsCut2)
print('mean =', mean_clip2)
print('stdev = ', stdev_clip2)
print(' ')
# Fit Gaussian
nbins = 100
xarray,yarray,mode,cmin,cmax =\
     gaussFit(points_clip2,nbins,mean_clip2,stdev_clip2)
# Plot Figure
plt.figure()
plt.title(r'2 Iterations of 3 Sigma Clip')
plt.xlabel(r'Photometric Redshift - Spectroscopic Redshift')
plt.ylabel(r'Number of Galaxies')
plt.hist(points_clip2,range=[cmin,cmax], bins=nbins);
plt.plot(xarray,yarray,color="red",linewidth=1.0,label='Gaussian Fit')
plt.axvline(x=mean_clip2,linewidth=1.0,color="yellow",label='mean')
plt.legend()
text_stats = r'''\noindent $\mu = {:.3f}$ \\'''.format(mean_clip2) + \
             r'''$\sigma = {:.3f}$ \\'''.format(stdev_clip2)
text_clip = r'''Total Clippied = {}'''.format(totclip)
#plt.text(0.4,150,text_stats,{'fontsize':20})
#plt.text(-0.6,190,text_clip)
if save:
    plt.savefig('SpecPhoto_clip2_{}.png'.format(tag),format='png',dpi=1000,bbox_inches='tight')

### Sigma Clipping 3 ###
points_clip3,pointsCut3,mean_clip3,stdev_clip3 = \
    sigma3clip(points_clip2,mean_clip2,stdev_clip2)
totclip = totclip + pointsCut3
print('3sigma clipped, iter3:',pointsCut3)
print('mean =', mean_clip3)
print('stdev = ', stdev_clip3)
# Fit Gaussian
xarray,yarray,mode,cmin,cmax =\
    gaussFit(points_clip3,nbins,mean_clip3,stdev_clip3)
# Plot Figure
plt.figure()
plt.title(r'3 Iterations of 3 Sigma Clip')
plt.xlabel(r'Photometric Redshift - Spectroscopic Redshift')
plt.ylabel(r'Number of Galaxies')
plt.hist(points_clip3,range=[cmin,cmax], bins=nbins);
plt.plot(xarray,yarray,color="red",linewidth=1.0,label='Gaussian Fit')
plt.axvline(x=mean_clip3,linewidth=1.0,color="yellow",label='mean')
plt.legend()
text_stats = r'''\noindent $\mu = {:.3f}$ \\'''.format(mean_clip3) + \
       r'''$\sigma = {:.3f}$ \\'''.format(stdev_clip3)
text_clip = r'''Total Clippied = {}'''.format(totclip)
#plt.text(0.2,100,text_stats,{'fontsize':20})
#plt.text(-0.4,120,text_clip)
if save:
    plt.savefig('SpecPhoto_clip3_{}.png'.format(tag),format='png',dpi=1000,bbox_inches='tight')

if show:
    plt.show()
plt.clf



