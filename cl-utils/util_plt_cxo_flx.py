'''
2020-10-04, Dennis Alp, dalp@kth.se

Plots aperture photometry from the CXO CSC 2.0 (products from CSCView).
'''

from __future__ import division, print_function
import os
from pdb import set_trace as db
import sys
import time
from glob import glob
from datetime import date
import subprocess

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import scipy.stats as sts
from scipy.ndimage import gaussian_filter
from astropy.coordinates import search_around_sky
from astropy.stats import LombScargle
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.time import Time

ff = glob(sys.argv[1])

for f in ff:
    print(f)
    f = fits.open(f)[4].data
    x = f['intensity']
    y = f['MPDF']
    plt.plot(x, y)

plt.show()
