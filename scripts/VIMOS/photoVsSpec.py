"""
Emily Biermann
VIMOS Analysis
Photo vs Spec z
6/5/19
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
cat3 = False

# Data Directories

pwd  = "/home/ebiermann/cat/VIMOS/"
if cat1:
    fname_match = 'MACS0329_match5b.fits'
if cat2:
    fname_match = 'MACS1347_match5b.fits'
if cat3:
    fname_match = 'MACS2211_match5b.fits'

specZ_idx = int(10)  # z_found
photZ_idx = int(13) # BPZ_Z_B

# Save Figures?
save = False
if cat1:
    figpwd = '/home/ebiermann/Cosmo-Research/figures/VIMOS/MACS0329/'
    tag = 'MACS0329' # tag for figure names
if cat2:
    figpwd = '/home/ebiermann/Cosmo-Research/figures/VIMOS/MACS1347/'
    tag = 'MACS1347' # tag for figure names
if cat3:
    figpwd = '/home/ebiermann/Cosmo-Research/figures/VIMOS/MACS2211/'
    tag = 'MACS2211' # tag for figure names

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

print 'Catalog: ' + fname_match
print 'Located in ' + pwd
print
print 'Notes: '

matchCat = fits.open(pwd + fname_match)
data = matchCat[1].data
matchCat.close()

specZ = []
photoZ = []

for i in range(0,len(data)):
    specZ.append(data[i][specZ_idx])
    photoZ.append(data[i][photZ_idx])

print 'STATISTICS'
print 
print 'Number of data points = {}'.format(len(specZ))
print

# Scatter plot with 1-1 line
x=np.linspace(0,np.amax(specZ))
y=x

plt.figure()
plt.title(r'Redshift Comparison for {}'.format(tag))
plt.xlabel(r'SpecZ from VIMOS')
plt.ylabel(r'PhotoZ from WtG')
#plt.errorbar(specZ, photoZ, xerr=specZ_err, yerr=photoZ_err, fmt='x')
plt.errorbar(specZ, photoZ,fmt='.')
plt.plot(x,y,linestyle='--')
if save:
    plt.savefig(figpwd+'SpecPhoto_scatter_{}.png'.format(tag),format='png',dpi=1000,bbox_inches='tight')

# Histogram
nbins = 50

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
print 'Unclipped'
print 'mean = {}'.format(mean)
print 'stdev = {}'.format(stdev)
print

# Plot Figure
plt.figure()
plt.xlim(-1.2,1.2)
plt.title(r'Redshift Comparison for {}'.format(tag))
plt.xlabel(r'Photometric Redshift - Spectroscopic Redshift')
plt.ylabel(r'Number of Galaxies')
plt.hist(points,range=[cmin,cmax], bins=nbins);
#plt.plot(xarray,yarray,color="red",linewidth=1.0,label='Gaussian Fit')
#plt.axvline(x=mean,linewidth=1.0,color="yellow",label='mean')
#plt.legend()
#plt.text(1.3,325,text_stats,{'fontsize':20})
if save:
    plt.savefig(figpwd+'SpecPhoto_noclip_{}.png'.format(tag),format='png',dpi=1000,bbox_inches='tight')

## 1+z correction and stats
cut = 0.1

points_all = points
points = []
cutpoints=0
i=0
for pZ in photoZ:
    point = (pZ - specZ[i])/(1+specZ[i])
    if abs(point) <= 0.1:
        points.append(point)
    else:
        cutpoints = cutpoints+1
    i=i+1

mean = np.average(points)
stdev = np.std(points)
xarray,yarray,mode,cmin_new,cmax_new =\
     gaussFit(points,nbins,mean,stdev)
print '1+z Corrected, DeltaZ > 0.1'
print 'mean = {}'.format(mean)
print 'stdev = {}'.format(stdev)
print 'fraction cut = {}'.format(cutpoints/float(len(points_all)))
print 
# Plot Figure
plt.figure()
#plt.xlim(-1.2,1.2)
plt.title(r'Redshift Comparison for {}'.format(tag))
plt.xlabel(r'$(\textrm{z}_\textrm{phot} - \textrm{z}_\textrm{spec})/(1+\textrm{z}_\textrm{spec})$')
plt.ylabel(r'Number of Galaxies')
plt.hist(points,range=[cmin_new,cmax_new], bins=nbins,label=r'$\Delta z$ Distribution');
plt.plot(xarray,yarray,color="red",linewidth=1.0,label='Gaussian Fit')
#plt.axvline(x=mean,linewidth=1.0,color="yellow",label='mean')
plt.legend()
text_stats = r'''\noindent $\mu = {:.3f}$ \\'''.format(mean) + \
             r'''$\sigma = {:.3f}$ \\'''.format(stdev)
#plt.text(0.05,38.0,text_stats,{'fontsize':15})
if save:
    plt.savefig(figpwd+'SpecPhoto_zCorrected_{}.png'.format(tag),format='png',dpi=1000,bbox_inches='tight')

'''
### Sigma Clipping 1 ###
points_clip1,pointsCut1,mean_clip1,stdev_clip1 = sigma3clip(points,mean,stdev)
totclip = pointsCut1
print '3sigma clipped, iter1: {}'.format(pointsCut1)
print 'mean = {}'.format(mean_clip1)
print 'stdev = {}'.format(stdev_clip1)
print 'fraction clipped = {}'.format(float(totclip)/float(len(specZ)))
print 
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
#plt.text(0.4,175,text_stats,{'fontsize':20})
#plt.text(-0.6,225,text_clip)
if save:
    plt.savefig(figpwd+'SpecPhoto_clip1_{}.png'.format(tag),format='png',dpi=1000,bbox_inches='tight')

### Sigma Clipping 2 ###
points_clip2,pointsCut2,mean_clip2,stdev_clip2 = \
    sigma3clip(points_clip1,mean_clip1,stdev_clip1)
totclip = totclip + pointsCut2
print '3sigma clipped, iter2: {}'.format(pointsCut2)
print 'mean = {}'.format(mean_clip2)
print 'stdev = {}'.format(stdev_clip2)
print 'fraction clipped = {}'.format(float(totclip)/float(len(specZ)))
print 
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
#plt.text(0.4,150,text_stats,{'fontsize':20})
#plt.text(-0.6,190,text_clip)
if save:
    plt.savefig(figpwd+'SpecPhoto_clip2_{}.png'.format(tag),format='png',dpi=1000,bbox_inches='tight')

### Sigma Clipping 3 ###
points_clip3,pointsCut3,mean_clip3,stdev_clip3 = \
    sigma3clip(points_clip2,mean_clip2,stdev_clip2)
totclip = totclip + pointsCut3
print '3sigma clipped, iter3: {}'.format(pointsCut3)
print 'mean = {}'.format(mean_clip3)
print 'stdev = {}'.format(stdev_clip3)
print 'fraction clipped = {}'.format(float(totclip)/float(len(specZ)))
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
#plt.text(0.2,100,text_stats,{'fontsize':20})
#plt.text(-0.4,120,text_clip)
if save:
    plt.savefig(figpwd+'SpecPhoto_clip3_{}.png'.format(tag),format='png',dpi=1000,bbox_inches='tight')
'''
if show:
    plt.show()

plt.clf()



