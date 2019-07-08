"""
Emily Biermann
Convert cat to fits_ldac format for scamp
6/18/19
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import stats
from scipy.stats import norm

pwd = '/home/ebiermann/cat/MACS0454_Astrometry/scamp/'
fname = 'crawford_spec_SCAMP_ref13.fits'


