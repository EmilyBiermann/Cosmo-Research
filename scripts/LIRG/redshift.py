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

#-------------------------------------------------------------------------------

pwd = '/home/ebiermann/cat/LIRG/'
fname = 'ms0451_crawSpec_match.fts'

# Save Figures?
save = False
if save:
    figpwd = '/home/ebiermann/Cosmo-Research/figures/LIRG/'
    tag = 'all'

# Show Figures?
show = True

#-------------------------------------------------------------------------------


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

#-------------------------------------------------------------------------------

print 'Catalog: ' + fname
print 'Located in ' + pwd
print
print 'Notes: '


matchCat = fits.open(pwd + fname)
data = matchCat[1].data
#cols = matchCat[1].columns
matchCat.close()

Z = data['Z']

nbins=25
mean = np.average(Z)
stdev = np.std(Z)
xarray,yarray,mode,cmin,cmax =\
     gaussFit(Z,nbins,mean,stdev)
print 'Unclipped'
print 'mean = {}'.format(mean)
print 'stdev = {}'.format(stdev)
print

# Plot Figure
plt.figure()
#plt.xlim(-1.2,1.2)
plt.title(r'FIR Redshift')
plt.xlabel(r'SpecZ')
plt.ylabel(r'Number of Galaxies')
plt.hist(Z,range=[cmin,cmax], bins=nbins);
#plt.plot(xarray,yarray,color="red",linewidth=1.0,label='Gaussian Fit')
#plt.axvline(x=mean,linewidth=1.0,color="yellow",label='mean')
#plt.legend()
#plt.text(1.3,325,text_stats,{'fontsize':20})
if save:
    plt.savefig(figpwd+'SpecZ_hist_{}.png'.format(tag),format='png',dpi=1000,bbox_inches='tight')

if show:
    plt.show()
plt.clf()

