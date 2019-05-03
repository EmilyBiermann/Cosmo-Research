"""
Emily Biermann
QQ Plot
4/23/19
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

# Path to catalogs
pwd = "/home/ebiermann/cat/matchCat/"
fname_match = 'WtG_SpecCrawNonzeroZ_13as.fits'

# Save Figures?
save = True
if save:
    figpwd = '/home/ebiermann/Cosmo-Research/figures/matchCat/all/QQ/'
    tag = 'MACS0454_allData' # tag for figure name

# Show Figures?
show = True

#-------------------------------------------------------------------------------

# Definitions

# QuanVal determines quantile of value [a] given numberline [x], P(x) [p]
# and step size [step]
def QuanVal(x,p,a,step):
    sum = 0.0
    for i in range(0,len(x)):
        sum = sum + step*p[i]
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

step = 0.01

# plot n random plots
nplots = 10
plotNums=np.random.randint(0,len(data),size=nplots)
for i in range(0,len(data)):
    specZ = data[i][776] # Z, 777
    pdz = data[i][662]   # pdz, 663
    z = np.arange(0.01,4.01,step)
    q = QuanVal(z,pdz,specZ,step)
    quant.append(q)
    # plot P(z) for random objects
    if np.isin(i,plotNums):
        SeqNr = data[i][2] # SeqNr_1, 3
        Rmag = data[i][650] # MAG_AUTO-SUBARU-COADD-1-W-C-RC, 651
        plt.figure()
        plt.title('P(z) for {}, Rmag = {:.2f}'.format(SeqNr,Rmag))
        plt.xlabel('z')
        plt.ylabel('P(z)')
        plt.plot(z,pdz,label=r'P(z) Distribution')
        plt.axvline(x=specZ,color='orange',label=r'Spectroscopic Redshift')
        plt.legend()
        if save:
            plt.savefig(figpwd+'pzPlot_{}_{}.png'.format(SeqNr,tag),\
            format='png',dpi=1000,bbox_inches='tight')
    else:
        continue


# QQ plot
nbins = 50

plt.figure()
plt.title(r'MACS0454')
plt.ylabel('n(Q)')
plt.xlabel('Q')
plt.hist(quant,bins=nbins)
    
if save:
    plt.savefig(figpwd+'QQplot_{}.png'.format(tag),format='png',dpi=1000,bbox_inches='tight')
if show:
    plt.show()
plt.clf()
