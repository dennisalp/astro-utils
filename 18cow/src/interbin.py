'''
2020-06-02, Dennis Alp, dalp@kth.se

Simulate the null distribution for Fourier interbins.
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
from astropy.io import fits
from astropy.time import Time
import scipy.stats as sts

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




################################################################
nn = 100000000
yy = sts.norm(0,1).rvs(nn)
ft = np.fft.rfft(yy)
po = np.abs(ft)**2/nn
ip = 1/np.sqrt(2)*(ft[1:-1]-ft[2:])
# Yes, of course it's 1/sqrt(2)!
ip = np.abs(ip)**2/nn
nu = np.fft.rfftfreq(nn, 1)

xx=np.arange(100)
plt.plot(xx, np.exp(-xx))
plt.hist(po, xx, normed=True, histtype='step')
plt.hist(ip, xx, normed=True, histtype='step')
plt.gca().set_yscale('log')
plt.show()
