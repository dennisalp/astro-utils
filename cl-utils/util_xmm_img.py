#!/Users/silver/anaconda3/envs/astroconda/bin/python

'''
2020-10-15, Dennis Alp, dalp@kth.se

Put this script in XMM PPS folders to create images.

It takes event files from EPIC and makes images.

Useful for cuts in energy or time (images not integrated over the entire observation and energies).
'''

from __future__ import division, print_function
import os
from pdb import set_trace as db
import sys
import time
from glob import glob
from datetime import date

import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from scipy.ndimage import gaussian_filter
from astropy import units as u
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units

def get_wcs(ff):
    dd = fits.open(ff)
    wcs = WCS()
    wcs.wcs.crpix = [dd[0].header['REFXCRPX'], dd[0].header['REFYCRPX']]
    wcs.wcs.cdelt = [dd[0].header['REFXCDLT'], dd[0].header['REFYCDLT']]
    wcs.wcs.crval = [dd[0].header['REFXCRVL'], dd[0].header['REFYCRVL']]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return wcs


################################################################
files = ['P0802860201M1S002MIEVLI0000.FTZ', 'P0802860201M2S003MIEVLI0000.FTZ', 'P0802860201PNS001PIEVLI0000.FTZ']
cams = ['M1', 'M2', 'PN']
img_size = 1024 # pixel
XMMEA_EP = 65584
XMMEA_EM = 65000


################################################################
xx = np.empty(0)
yy = np.empty(0)

for ii, ff in enumerate(files):
    cam = cams[ii]
    dd = fits.open(ff)[1].data

    # Concerning flags
    # https://xmm-tools.cosmos.esa.int/external/sas/current/doc/eimageget/node4.html
    # http://xmm.esac.esa.int/xmmhelp/EPICpn?id=8560;expression=xmmea;user=guest
    if cam == 'PN':
        good = dd['FLAG'] <= XMMEA_EP
        good = good & (dd['PATTERN'] <= 4)
    elif cam == 'M1' or cam == 'M2':
        good = dd['FLAG'] <= XMMEA_EM
        good = good & (dd['PATTERN'] <= 12)
    else:
        print('ERROR Unknown instrument:', cam, ii)
        sys.exit(1)
        
    good = good & (dd['PI'] > 300.) & (dd['PI'] < 9500.)
    good = good & (dd['X'] > 0.) # XMM puts some invalid photons at (-99999999.0, -99999999.0)
    good = good & (dd['Y'] > 0.)
    xx = np.concatenate((xx, dd['X'][good]))
    yy = np.concatenate((yy, dd['Y'][good]))

wcs = get_wcs(ff)
img, xbin, ybin = np.histogram2d(xx, yy, img_size)
img = np.log10(img+0.1).T

plt.imshow(img, origin='lower')
plt.figure()
plt.imshow(gaussian_filter(img, 1), origin='lower')
plt.show()
db()
