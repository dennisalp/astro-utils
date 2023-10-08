'''
2020-10-03, Dennis Alp, dalp@kth.se

HLX catalog of barrows19. Sky images from XMM, CXO, XRT,
PanSTARRS, DSS, and HST/ACS.
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
def get_dat(ff):
    ff = open(ff, 'r')
    RA = []
    DE = []
    for ll in ff.readlines():
        if ll[0] == '#':
            continue
        tmp = ll.split()[5]
        ra = ' '.join([tmp[1:3], tmp[3:5], tmp[5:9]])
        de = ' '.join([tmp[9:12], tmp[12:14], tmp[14:16]])
        RA.append(ra)
        DE.append(de)

    src = SkyCoord(RA, DE, unit=(u.hourangle, u.deg))
    return src

################################################################
src = get_dat('/Users/silver/Box Sync/sci/lib/b/barrows19/table1.txt')
cwd = '/Users/silver/box/sci/doc/proposals/20_xmm_ulx/out/barrows19/'
os.chdir(cwd)

################################################################
for ii, ss in enumerate(src):
    mk_plt(ss)
    plt.savefig('{0:d}.pdf'.format(ii), bbox_inches='tight', pad_inches=0.1, dpi=300)
    # plt.show()

db()
