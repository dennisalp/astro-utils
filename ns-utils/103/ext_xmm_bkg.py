'''
2020-05-27, Dennis Alp, dalp@kth.se

Prepare the reduced and barycentered HRC-S event lists into pure timestamps.
'''

from __future__ import division, print_function
import os
from pdb import set_trace as db
import sys
from glob import glob
import time
from datetime import date

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.time import Time




################################################################
def mk_inf(oi, t0, tt, tb, dt):
    ff  = oi
    obj = '1E 161348-5055'
    ra  = '16:17:36.23'
    de  = '+51:02:24.6s'
    mjd = Time(t0, format='cxcsec').mjd
    nn  = tt.size
    ww  = dt
    obs = 'Dennis Alp'

    inf = open(out + '{0:s}.inf'.format(ff), 'w')
    buf = """ Data file name without suffix          =  {0:s}
 Telescope used                         =  XMM
 Instrument used                        =  EPIC/pn
 Object being observed                  =  {1:s}
 J2000 Right Ascension (hh:mm:ss.ssss)  =  {2:s}
 J2000 Declination     (dd:mm:ss.ssss)  =  {3:s}
 Data observed by                       =  XMM
 Epoch of observation (MJD)             =  {4:f}
 Barycentered?           (1=yes, 0=no)  =  1
 Number of bins in the time series      =  {5:d}
 Width of each time series bin (sec)    =  {6:f}
 Any breaks in the data? (1=yes, 0=no)  =  0
 Type of observation (EM band)          =  X-ray
 Field-of-view diameter (arcsec)        =  2.0
 Central energy (kev)                   =  2.0
 Energy bandpass (kev)                  =  4.0
 Data analyzed by                       =  {7:s}
 Any additional notes:
    None
""".format(ff, obj, ra, de, mjd, nn, ww, obs)

    inf.write(buf)
    inf.close()



def do_stuff(tt):
    np.savetxt(out + oi + '.txt', tt)

    t0 = tt[0]
    tt = tt-t0
    tb = np.arange(tt[0], tt[-1]+dt, dt)
    tt, _ = np.histogram(tt, tb)
    mk_inf(oi, t0, tt, tb, dt)
#    tt.astype(np.float64).tofile(out + oi + '.bin')


################################################################
files = sorted(glob('/Users/silver/dat/xmm/103/0743750201_repro/epn_cl_bkg_bc.evt'))
out = '/Users/silver/box/phd/pro/nss/103/dat/'
dt = 1.e-3

for ii, ff in enumerate(files):
    oi = ff.split('/')[6].replace('repro', 'bkg')
    dat = fits.open(ff)[1].data
    tt = dat['time']
    # REMOVE INSTRUMENTAL ARTIFACTS, CRITICALLY IMPORTANT
    tt = tt + np.random.uniform(0, 0.005671, tt.size)
    do_stuff(tt)

db()
