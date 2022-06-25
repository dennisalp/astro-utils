'''
0000-00-00, Dennis Alp, dalp@kth.se

Description of this file.
'''

from __future__ import division, print_function
import os
from pdb import set_trace as db
import sys
from glob import glob
import time
from datetime import date
from tqdm import tqdm

import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units
from astropy.time import Time

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

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
mp = 1.67262192369e-24 # g

sna = SkyCoord(83.86661853961974*units.degree, -69.26975365456684*units.degree, frame='icrs') # SN 1987A, my note, 5:35:27.9884, -69:16:11.1132
cas = SkyCoord('23h23m27.943s', '+58d48m42.51s', frame='icrs') # Cas A CCO, fesen06, (350.86642917, 58.81180833)
cow = SkyCoord('16h16m0.22418s', '+22d16m04.8903s', frame='icrs') # AT 2018cow, bietenholz20, (244.00093408, 22.26802508)
snj = SkyCoord('09h55m42.137s', '+69d40m25.40s', frame='fk5') # SN 2014J, kelly14
rcw = SkyCoord('16h17m36.23s', '+51d02m24.6s', frame='fk5') # 1E 161348-5055, RCW 103, de_luca08




################################################################
tt = np.empty(0)

print('Reading')
for ff in sys.argv[1:]:
    tt = np.concatenate((tt, fits.open(ff)[1].data['TIME']))

tt = np.sort(tt)
tt = tt[:2433135]
nn = 27
dt = 0.0005
nu = np.fft.rfftfreq(2**nn, dt)
tt = tt-tt[0]

print('Binning')
if True:
    tmp = np.zeros(int(np.ceil(tt[-1]/dt)))
    tmp[(tt//dt).astype(np.int)] += 1 # assumes only 1 photon per bin
else:
    tmp, _ = np.histogram(tt, np.arange(tt[0], tt[-1], dt))

lc = np.ones(2**nn)*tmp.mean()
lc[:np.min((tmp.size, lc.size))] = tmp[:np.min((tmp.size, lc.size))]

print('FFTing')
ft = np.fft.rfft(lc)
po = np.abs(ft)**2/tt.size

ii = po > 18
plt.loglog(nu[ii], po[ii])
plt.show()
db()
