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

# Data Directories
# Topcat match MUST have WtG as first catalog!!

pwd  = "/home/ebiermann/cat/mag_matchCat/"
if cat1:
    fname_match = 'match4_1as1mag_3as1mag.fits'
if cat2:
    fname_match = 'match3_3as1mag_rem.fits'

# Save Figures?
save = True
if save:
    if cat1:
        figpwd = '/home/ebiermann/Cosmo-Research/figures/mag_matchCat/all/QQ/'
        tag = 'all' # tag for figure names
    if cat2:
        figpwd = '/home/ebiermann/Cosmo-Research/figures/mag_matchCat/match3_rem/QQ/'
        tag = 'match3rem'

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
        #plt.title('P(z) for {}, Rmag = {:.2f}'.format(SeqNr,Rmag))
        plt.title('P(z) Distribution, Rmag = {:.2f}'.format(Rmag))
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
    
    


# PIT/QQ residuals Plot
nbins = 75
Qdata = np.sort(quant)
Qtheory = np.linspace(0.0,1.0,len(quant))
delQ = Qdata - Qtheory

plt.figure()
gridspec.GridSpec(3,2)

plt.subplot2grid((3,2), (0,0), colspan=2, rowspan=2)
plt.title('MACS0454')
#plt.xlabel('PIT Value')
plt.ylabel('Number of Galaxies')
plt.hist(quant,bins=nbins)

plt.subplot2grid((3,2), (2,0), colspan=2, rowspan=1)
plt.xlabel(r'Quantile')
plt.ylabel(r'$\Delta_{\textrm{Q}}$')
plt.hlines(0.0,0.0,1.0,linestyles='--')
plt.plot(Qtheory,delQ)

plt.tight_layout()
if save:
    plt.savefig(figpwd+'QQplot_{}.png'.format(tag),\
    format='png',dpi=1000,bbox_inches='tight')


if show:
    plt.show()
plt.clf()
