'''
2020-10-03, Dennis Alp, dalp@kth.se

Find all sources in 4XMM-DR9 (webb20) and match with sources from
earnshaw19b. This is just a quick and dirty way to find all 4XMM-DR9
source in nearby galaxies. Append some of the following galaxies (not
all) that were excluded by earnshaw19b: Andromeda, the Large and Small
Magellanic Clouds, M33, M54, M81, M101, NGC 55, NGC 253, NGC 5128,
Draco Dwarf, Sculp- tor Dwarf, and Sextans Dwarf Spheroidal.

Ignore MCs. M31 covered by stiele11 and nothing interesting in M33
based on long10.

Produce light XMM, CXO, and XRT curves and sky images from XMM, CXO,
XRT, PanSTARRS, DSS, and HST/ACS for all sources.

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
from xmm_lc import mk_plt

################################################################
def get_xmm(ff):
    dat = fits.open(ff)[1].data
    flx = dat['EP_8_FLUX']
    idx = np.where(flx > xmm_flx_lim)[0]
    flx = flx[idx]
    src = dat['SRCID'][idx]
    ra = dat['SC_RA'][idx]
    de = dat['SC_DEC'][idx]
    coo = SkyCoord(ra, de, unit=(u.deg, u.deg))
    return src, flx, coo

def get_ulx(ff):
    dat = fits.open(ff)[1].data
    ulx = dat['SRCID']
    ulx, idx = np.unique(ulx, return_index=True)
    dis = dat['DISTANCE_BEST'][idx].astype(np.float64)
    ra = dat['SC_RA'][idx]
    de = dat['SC_DEC'][idx]
    coo = SkyCoord(ra, de, unit=(u.deg, u.deg))
    return ulx, dis, coo

################################################################
xmm_flx_lim = 3.e-13 # Keep sources with a flux higher than this in at least 1 detection

# Large galaxies excluded by earnshaw19b
# M81, M101, NGC 55, NGC 253, NGC 5128 (Cen A)
ra = ['09 55 33.173', '14 03 12.583', '00 14 53.602', '00 47 33.134', '13 25 27.615']
de = ['+69 03 55.06', '+54 20 55.50', '-39 11 47.86', '-25 17 19.68', '-43 01 08.80']
close_galaxies = SkyCoord(ra, de, unit=(u.hourangle, u.deg))

xmm_src, flx, xmm_coo = get_xmm('/Users/silver/box/sci/lib/w/webb20/4XMM_DR9cat_v1.0.fits')
ulx, dis, ulx_coo = get_ulx('/Users/silver/box/sci/lib/e/earnshaw19b/catalog.fits')
cwd = '/Users/silver/box/sci/doc/proposals/20_xmm_ulx/out/webb20/'
os.chdir(cwd)


################################################################
# Find coordinates from 4XMM-DR9
ulx_idx, xmm_idx, _, _ = search_around_sky(ulx_coo, xmm_coo, 10*u.arcmin)
_, xmm_idx_2, _, _ = search_around_sky(close_galaxies, xmm_coo, 45*u.arcmin)
xmm_idx = np.concatenate((xmm_idx, xmm_idx_2))

xmm_src = xmm_src[xmm_idx]
xmm_coo = xmm_coo[xmm_idx]
_, idx = np.unique(xmm_src, return_index=True)
src = xmm_coo[idx]

# Make plots for all selected coordinates
for ii, ss in enumerate(src):
    print(ii, src.size)
    if ii < int(sys.argv[1]):
        continue
    
    mk_plt(ss)
    plt.savefig('{0:d}.pdf'.format(ii), bbox_inches='tight', pad_inches=0.1, dpi=300)
    # plt.show()
    plt.close()

# db()
