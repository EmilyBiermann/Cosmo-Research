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

cat1 = False
cat2 = True

# Data Directories
# Topcat match MUST have WtG as first catalog!!

pwd  = "/home/ebiermann/cat/mag_matchCat/"
if cat1:
    fname_match = 'match4_1as_3as1mag.fits'
if cat2:
    fname_match = 'match3_3as1mag_rem.fits'

# Save Figures?
save = False
if save:
    if cat1:
        pwdfig = '/home/ebiermann/Cosmo-Research/figures/mag_matchCat/all/mag/'
        tag = 'all' # tag for figure names
    if cat2:
        pwdfig = '/home/ebiermann/Cosmo-Research/figures/mag_matchCat/match3_rem/mag/'
        tag = 'match3rem'

# Show Figures?
show = True

#-------------------------------------------------------------------------------

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

matchCat = fits.open(pwd + fname_match)
data = matchCat[1].data

WtG_mag = []
CRAW_mag = []

for i in range(0,len(data)):
    WtG_mag.append(data[i][650])  # MAG_AUTO-SUBARU-COADD-1-W-C-RC
    CRAW_mag.append(data[i][674]) # MAG, 675


# 1-1 line
x=np.linspace(0,np.amax(CRAW_mag))
y=x

plt.figure()
plt.title(r'Magnitude Comparison')
plt.xlabel(r'Mag from Crawford')
plt.ylabel(r'Mag from WtG')
plt.errorbar(CRAW_mag, WtG_mag, fmt='.')
plt.plot(x,y,linestyle='--')
if save:
    plt.savefig(pwdfig+'MagMatch_line_{}.png'.format(tag),format='png',dpi=1000,bbox_inches='tight')

# Histogram
nbins = 100

# Unclipped Data
points = []
i=0
for m in WtG_mag:
    points.append(m - CRAW_mag[i])
    i=i+1
mean = np.average(points)
stdev = np.std(points)
print('Unclipped')
print('mean =', mean)
print('stdev = ', stdev)
print(' ')
xarray,yarray,mode,cmin,cmax =\
     gaussFit(points,nbins,mean,stdev)
# Plot Figure
plt.figure()
plt.title(r'Unclipped Data')
plt.xlabel(r'WtG Mag - Crawford Mag')
plt.ylabel(r'Number of Galaxies')
plt.hist(points,range=[cmin,cmax], bins=nbins);
plt.plot(xarray,yarray,color="red",linewidth=1.0,label='Gaussian Fit')
plt.axvline(x=mean,linewidth=1.0,color="yellow",label='mean')
plt.legend()
text_stats = r'''\noindent $\mu = {:.3f}$ \\'''.format(mean) + \
             r'''$\sigma = {:.3f}$ \\'''.format(stdev)
#plt.text(39,900,text_stats,{'fontsize':20})
if save:
    plt.savefig(pwdfig+'MagMatch_noclip_{}.png'.format(tag),format='png',dpi=1000,bbox_inches='tight')

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
plt.xlabel(r'WtG Mag - Crawford Mag')
plt.ylabel(r'Number of Galaxies')
plt.hist(points_clip1,range=[cmin,cmax], bins=nbins);
plt.plot(xarray,yarray,color="red",linewidth=1.0,label='Gaussian Fit')
plt.axvline(x=mean_clip1,linewidth=1.0,color="yellow",label='mean')
plt.legend()
text_stats = r'''\noindent $\mu = {:.3f}$ \\'''.format(mean_clip1) + \
             r'''$\sigma = {:.3f}$ \\'''.format(stdev_clip1)
text_clip = r'''Total Clippied = {}'''.format(totclip)
#plt.text(1.7,375,text_stats,{'fontsize':20})
#plt.text(-12.5,490,text_clip)
if save:
    plt.savefig(pwdfig+'MagMatch_clip1_{}.png'.format(tag),format='png',dpi=1000,bbox_inches='tight')

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
plt.xlabel(r'WtG Mag - Crawford Mag')
plt.ylabel(r'Number of Galaxies')
plt.hist(points_clip2,range=[cmin,cmax], bins=nbins);
plt.plot(xarray,yarray,color="red",linewidth=1.0,label='Gaussian Fit')
plt.axvline(x=mean_clip2,linewidth=1.0,color="yellow",label='mean')
plt.legend()
text_stats = r'''\noindent $\mu = {:.3f}$ \\'''.format(mean_clip2) + \
             r'''$\sigma = {:.3f}$ \\'''.format(stdev_clip2)
text_clip = r'''Total Clippied = {}'''.format(totclip)
#plt.text(2.6,400,text_stats,{'fontsize':20})
#plt.text(-7,490,text_clip)
if save:
    plt.savefig(pwdfig+'MagMatch_clip2_{}.png'.format(tag),format='png',dpi=1000,bbox_inches='tight')

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
plt.xlabel(r'WtG Mag - Crawford Mag')
plt.ylabel(r'Number of Galaxies')
plt.hist(points_clip3,range=[cmin,cmax], bins=nbins);
plt.plot(xarray,yarray,color="red",linewidth=1.0,label='Gaussian Fit')
plt.axvline(x=mean_clip3,linewidth=1.0,color="yellow",label='mean')
plt.legend()
text_stats = r'''\noindent $\mu = {:.3f}$ \\'''.format(mean_clip3) + \
       r'''$\sigma = {:.3f}$ \\'''.format(stdev_clip3)
text_clip = r'''Total Clippied = {}'''.format(totclip)
#plt.text(1.3,250,text_stats,{'fontsize':20})
#plt.text(-3,345,text_clip)
if save:
    plt.savefig(pwdfig+'MagMatch_clip3_{}.png'.format(tag),format='png',dpi=1000,bbox_inches='tight')

if show:
    plt.show()
plt.clf


