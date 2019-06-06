"""
Emily Biermann
QQ Plot
4/23/19
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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
cat3 = False

# Data Directories

pwd  = "/home/ebiermann/cat/VIMOS/"
if cat1:
    fname_match = 'MACS0329_match5b.fits'
if cat2:
    fname_match = 'MACS1347_match5b.fits'
if cat3:
    fname_match = 'MACS2211_match5b.fits'

specZ_idx = int(10) # z_found
photZ_idx = int(13) # BPZ_Z_B
pdz_idx = int(17)   # pdz

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

# Definitions

# QuanVal determines quantile of value [a] given numberline [x], P(x) [p]
# and step size [step]
def QuanVal(x,p,a,step):
    sum = 0.0
    for i in range(0,len(x)):
        sum = sum + step*p[i+1]
        if  a < x[i+1]:
            quant = sum
            break
        else:
            continue
    return(quant/step)

# Open catalog
cat = fits.open(pwd + fname_match)
data = cat[1].data
cat.close()

# Allocate arrays
quant = []
outlier = 0
step = 0.01

nbins=75
outLow = 1.0/float(nbins)
outHigh = 1.0 - outLow

# plot n random plots
nplots = 10
#plotNums=np.random.randint(0,len(data),size=nplots)
plotNums=(332,485,528)
for i in range(0,len(data)):
    specZ = data[i][specZ_idx] # z_found
    pdz = data[i][pdz_idx]   # pdz
    z = np.arange(0.01,4.01,step)
    q = QuanVal(z,pdz,specZ,step)
    if q <=outLow or q>=outHigh:
       outlier = outlier+1
    quant.append(q)
    '''
    # plot P(z) for random objects
    if np.isin(i,plotNums):
        SeqNr = data[i][2] # SeqNr_1, 3
        Rmag = data[i][650] # MAG_AUTO-SUBARU-COADD-1-W-C-RC, 651
        Z_ml = bpzData[SeqNr-1][6] # BPZ_Z_ML, 7
        Z_B = bpzData[SeqNr-1][1]
        plt.figure()
        plt.title('P(z) for {}, Rmag = {:.2f}'.format(SeqNr,Rmag))
        #plt.title('P(z) Distribution, Rmag = {:.2f}'.format(Rmag))
        plt.xlabel('z')
        plt.ylabel('P(z)')
        plt.xlim(0.0,1.5)
        plt.plot(z,pdz,label=r'P(z) Distribution')
        plt.axvline(x=Z_B,color='blue',linestyle='--',label=r'$Z_{B}$')  
        plt.axvline(x=Z_ml,color='green',linestyle=':',label=r'$Z_{ML}$')
        plt.axvline(x=specZ,color='orange',label=r'Spectroscopic Redshift')
        plt.legend()
        if save:
            plt.savefig(figpwd+'pzPlot_{}_{}.png'.format(SeqNr,tag),\
            format='png',dpi=1000,bbox_inches='tight')
    else:
        continue
    '''
    


# PIT/QQ residuals Plot

Qdata = np.sort(quant)
Qtheory = np.linspace(0.0,1.0,len(quant))
delQ = Qdata - Qtheory

nGal_th = float(len(quant))/nbins
fail = outlier - 2*nGal_th
print('Catastrophic Failures: {}'.format(fail))
print('Fraction of Cat. Fail: {}'.format(fail/float(len(quant))))

plt.figure()
gridspec.GridSpec(3,2)

plt.subplot2grid((3,2), (0,0), colspan=2, rowspan=2)
plt.title('{}'.format(tag))
#plt.xlabel('PIT Value')
plt.ylabel('Number of Galaxies')
plt.hist(quant,bins=nbins)
plt.hlines(nGal_th,0.0,1.0,colors='k',label=r'Even Distribution')

plt.subplot2grid((3,2), (2,0), colspan=2, rowspan=1)
plt.xlabel(r'Quantile')
plt.ylabel(r'$\Delta_{\textrm{Q}}$')
plt.ylim(-0.1,0.1)
plt.hlines(0.0,0.0,1.0,linestyles='--')
plt.plot(Qtheory,delQ)

plt.tight_layout()
if save:
    plt.savefig(figpwd+'QQplot_{}.png'.format(tag),\
    format='png',dpi=1000,bbox_inches='tight')


if show:
    plt.show()
plt.clf()
