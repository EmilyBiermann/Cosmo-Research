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

# Data Directories
# Topcat match MUST have WtG as first catalog!!

pwd  = "/home/ebiermann/cat/mag_matchCat/"
fname_match1 = 'match4_1as1mag_3as1mag.fits'
fname_match2 = 'match5_goodObjects.fits'

# Save Figures?
save = True
if save:
    figpwd = '/home/ebiermann/Cosmo-Research/figures/mag_matchCat/all/QQ/'
    tag = 'lensingAnalysis' # tag for figure names

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
cat = fits.open(pwd + fname_match1)
data = cat[1].data
cat.close()

cat = fits.open(pwd + fname_match2)
data2 = cat[1].data
cat.close()

# Allocate arrays
quant = []
quant2 = []

step = 0.01

for i in range(0,len(data)):
    specZ = data[i][776] # Z, 777
    pdz = data[i][662]   # pdz, 663
    z = np.arange(0.01,4.01,step)
    q = QuanVal(z,pdz,specZ,step)
    quant.append(q)

for i in range(0,len(data2)):
    specZ = data2[i][776] # Z, 777
    pdz = data2[i][662]   # pdz, 663
    z = np.arange(0.01,4.01,step)
    q = QuanVal(z,pdz,specZ,step)
    quant2.append(q)

# PIT/QQ residuals Plot
nbins = 75
nbins2 = 25

Qdata = np.sort(quant)
Qtheory = np.linspace(0.0,1.0,len(quant))
delQ = Qdata - Qtheory

Qdata2 = np.sort(quant2)
Qtheory2 = np.linspace(0.0,1.0,len(quant2))
delQ2 = Qdata2 - Qtheory2

plt.figure()
gridspec.GridSpec(3,2)

plt.subplot2grid((3,2), (0,0), colspan=2, rowspan=2)
plt.title('MACS0454')
#plt.xlabel('PIT Value')
plt.ylabel('Number of Galaxies')
plt.hist(quant,bins=nbins,label='All Objects')
plt.hist(quant2,bins=nbins,label='Lensing Analysis')

plt.legend()

plt.subplot2grid((3,2), (2,0), colspan=2, rowspan=1)
plt.xlabel(r'Quantile')
plt.ylabel(r'$\Delta_{\textrm{Q}}$')
plt.hlines(0.0,0.0,1.0,linestyles='--')
plt.plot(Qtheory,delQ,label='All Objects')
plt.plot(Qtheory2,delQ2,label='Lensing Analysis')

plt.tight_layout()
if save:
    plt.savefig(figpwd+'QQplot_overlay_{}.png'.format(tag),\
    format='png',dpi=1000,bbox_inches='tight')


if show:
    plt.show()
plt.clf()
