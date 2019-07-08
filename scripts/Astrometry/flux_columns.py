"""
Emily Biermann
Add flux/err columns from mag to CrawSpec
6/18/19
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import stats
from scipy.stats import norm

magZpt = 27.0

pwd = '/home/ebiermann/cat/MACS0454_Astrometry/scamp/'
fname = 'crawford_spec_SCAMP.fits'

hdul = fits.open(pwd+fname)
data = hdul[1].data

cols=hdul[1].columns
names = cols.names[20:113]

cols_new = []
cols_new.append(cols['ID'])

for name in names:
    if 'err' in name:
        sigMag = data[name]
        sigF = []
        for i in range(0,len(sigMag)):
            sigF.append(np.log(10.0)/2.5*np.abs(F[i])*sigMag[i])
        cols_new.append(fits.Column(name=name+'_Flux', format='D', array=sigF))
    else:
        mag = data[name]
        F = []
        for i in range(0,len(mag)):
            F.append(10.0**(-(mag[i]-magZpt)/2.5))
        cols_new.append(fits.Column(name=name+'_Flux',format='D',array=F))

t = fits.BinTableHDU.from_columns(cols_new)
t.writeto(pwd+'crawSpec_Flux.fits')
        
hdul.close()
