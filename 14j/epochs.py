'''
2019-06-26, Dennis Alp, dalp@kth.se

Just computes the epoch of NuSTAR observations of SN 2014J.
'''

from __future__ import division, print_function
import os
from pdb import set_trace as db
import sys
from glob import glob
from datetime import datetime as dt
from datetime import timedelta as td

import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units

# Constants, cgs
cc = 2.99792458e10 # cm s-1
GG = 6.67259e-8 # cm3 g-1 s-2
hh = 6.6260755e-27 # erg s
DD = 51.2 # kpc
pc = 3.086e18 # cm
kpc = 3.086e21 # cm
mpc = 3.086e24 # cm
kev2erg = 1.60218e-9 # erg keV-1
Msun = 1.989e33 # g
Lsun = 3.828e33 # erg s-1
Rsun = 6.957e10 # cm
Tsun = 5772 # K
uu = 1.660539040e-24 # g
SBc = 5.670367e-5 # erg cm-2 K-4 s-1
kB = 1.38064852e-16 # erg K-1

snj = dt(2014, 1, 14, 18)




################################################################
# Manual version
oid = [80002092002, 80002092004, 80002092006, 80002092007, 80002092008, 80002092009, 80002092011]
telapse = np.array([1.225862500037551E+05, 1.705267500113845E+05, 5.794903199999332E+05, 5.620065250018984E+05, 6.129942499423027E+04, 2.128494249942303E+05, 2.009156250094771E+05])
exp = np.array([65818, 89966, 309525, 306240, 33648, 114894, 110925])

exp = np.array([
    [5.969902849137734E+04, 5.978341401583719E+04],
    [7.738947305070031E+04, 7.694133950998796E+04],
    [2.713249751612501E+05, 2.709355935283264E+05],
    [2.599899406944642E+05, 2.582891128849272E+05],
    [2.649675212473892E+04, 2.685222205733449E+04],
    [9.818484749761189E+04, 9.775469420900149E+04],
    [9.823970971892483E+04, 9.698085981570886E+04]])
exp = np.average(exp, axis=1)

epo = [dt(2014, 1, 23, 12, 36, 7) + td(seconds=0.5*telapse[0]), # THIS IS WRONG, NEED LIVETIME OR SOMETHING
       dt(2014, 1, 25, 19, 31, 7) + td(seconds=0.5*telapse[1]),
       dt(2014, 1, 28, 12, 11, 7) + td(seconds=0.5*telapse[2]),
       dt(2014, 2,  4,  5, 41, 7) + td(seconds=0.5*telapse[3]),
       dt(2014, 2, 10, 18, 26, 7) + td(seconds=0.5*telapse[4]),
       dt(2014, 2, 11, 12,  6, 7) + td(seconds=0.5*telapse[5]),
       dt(2014, 3,  3, 16, 56, 7) + td(seconds=0.5*telapse[6])]

for ii, ee in enumerate(epo):
    tt = ee - snj
    ep = ee.strftime("%Y-%m-%d %H:%M:%S.%f")[:-5]
    tp = tt.days+tt.seconds/86400
    wp = exp[ii]/exp.sum()
    print('{3:<11d} Epoch: {0:21} Time:{1:9.5f} Weight:{2:8.5f} Exposure:{4:8.3f}'.format(ep, tp, wp, oid[ii], exp[ii]/1e3))



################################################################
# The automatic below this line
files = [
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092002/nu80002092002A01_sr_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092002/nu80002092002A01_bk_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092002/nu80002092002B01_sr_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092002/nu80002092002B01_bk_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092004/nu80002092004A01_sr_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092004/nu80002092004A01_bk_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092004/nu80002092004B01_sr_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092004/nu80002092004B01_bk_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092006/nu80002092006A01_sr_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092006/nu80002092006A01_bk_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092006/nu80002092006B01_sr_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092006/nu80002092006B01_bk_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092007/nu80002092007A01_sr_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092007/nu80002092007A01_bk_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092007/nu80002092007B01_sr_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092007/nu80002092007B01_bk_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092008/nu80002092008A01_sr_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092008/nu80002092008A01_bk_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092008/nu80002092008B01_sr_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092008/nu80002092008B01_bk_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092009/nu80002092009A01_sr_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092009/nu80002092009A01_bk_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092009/nu80002092009B01_sr_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092009/nu80002092009B01_bk_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092011/nu80002092011A01_sr_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092011/nu80002092011A01_bk_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092011/nu80002092011B01_sr_grp.pha',
    '/Users/silver/dat/nus/14j/' + sys.argv[1] + '/80002092011/nu80002092011B01_bk_grp.pha']

exp = 0.
for ff in files:
    if 'bk' in ff: continue
    exp += fits.getval(ff,'EXPOSURE', 1)/2.

ex = np.zeros(2)
for ii, ff in enumerate(files):
    if np.mod(ii,2) == 0: srscal = fits.getval(ff,'BACKSCAL', 1)
    else: bkscal = fits.getval(ff,'BACKSCAL', 1)
    ex[np.mod(ii,4)//2] = fits.getval(ff,'EXPOSURE', 1)

    te = fits.getval(ff,'TELAPSE', 1)
    da = fits.getval(ff,'DATE-OBS', 1)
    ee = dt(int(da.split('T')[0].split('-')[0]),
            int(da.split('T')[0].split('-')[1]),
            int(da.split('T')[0].split('-')[2]),
            int(da.split('T')[1].split(':')[0]),
            int(da.split('T')[1].split(':')[1]),
            int(da.split('T')[1].split(':')[2])) + td(seconds=0.5*te)

    tt = ee - snj
    ep = ee.strftime("%Y-%m-%d %H:%M:%S.%f")[:-5]
    tp = tt.days+tt.seconds/86400
    wp = ex.mean()/exp

    if np.mod(ii,4) == 3:
        print('{3:<11d} Epoch: {0:21} Time:{1:9.5f} Weight:{2:8.5f} Exposure:{4:8.3f} Backfrac:{5:10.5f}'.format(ep, tp, wp, int(ff.split('/')[-1][2:13]), ex.mean()/1e3, srscal/bkscal))

db()
