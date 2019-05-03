# Emily Biermann
# 2/9/19
#
# Script to match photo-z (Von der Linden) to spec-z (Crawford) catalogs 
# for cluser MS0451

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
#from scipy import stats
#from scipy.stats import norm

pwd = "../cat/matchCat/"
photoCat = fits.open(pwd + "MACS0454-03.W-C-RC.cat")
specCat = fits.open(pwd + "w05.master.fits")
pData = photoCat[1].data
sData = specCat[1].data
photoCat.close()
specCat.close()

sRA = []
sDEC = []
sRmag = []
sRmag_err = []
sZ = []

pRA = []
pDEC = []

for i in range(0,len(sData)):
    # Note indicies are index - 1 in crawford
    ZTYPE = sData[i][115]
    if ZTYPE == 2 continue # If photoz, skip
    sRA.append(sData[i][3])
    sDEC.append(sData[i][4])
    sRmag.append(sData[i][59])     # RT
    sRmag_err.append(sData[i][60])
    sZ.append(sData[i][113])       # Best Redshift

1arcsec = 1.0/3600.0 # one arcsec in deg

for i in range(0,len(pData)):
    pRA = pData[i][13]    # ALPHA_J2000
    pDEC = pData[i][14]   # DELTA_J2000
    pRmag = pData[i][650] # MAG_AUTO_SUBARU_COADD-1-W-C-RC
    


