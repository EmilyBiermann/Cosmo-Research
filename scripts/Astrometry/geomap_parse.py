"""
Emily Biermann
geomap parsing
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

pwd  = "/home/ebiermann/cat/MACS0454_Astrometry/"
fnames = ['MACS0454_CrawSpec12_XY.txt','MACS0454_CrawSpec13_XY.txt',
          'MACS0454_CrawSpec17_XY.txt','MACS0454_CrawSpec18_XY.txt',
          'MACS0454_CrawSpec19_XY.txt']

'''
for fname in fnames:
    with open(pwd + fname,'rb') as f:
        clean_lines = (line.replace(b' ',b',') for line in f)
        data = np.genfromtxt(clean_lines, delimiter=',')
    np.savetxt(fname,data,delimiter=',')
'''
'''
for fname in fnames:
    data = np.loadtxt(pwd+fname,delimiter=',')
    np.savetxt(fname,data,delimiter=' ')
'''

pwd = '/home/ebiermann/cat/LIRG/'
fname = 'ms0451_crawSpec13_geomapDatabase.txt'
with open(pwd + fname,'rb') as f:
    clean_lines = (line.replace(b',',b' ') for line in f)
    data = np.genfromtxt(clean_lines, delimiter=' ')
    np.savetxt(fname,data,delimiter=' ')
