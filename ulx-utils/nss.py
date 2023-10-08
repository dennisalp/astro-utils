'''
2020-09-21, Dennis Alp, dalp@kth.se

Worksheet for pulsar luminosities and evolution.
'''

from __future__ import division, print_function
import os
from pdb import set_trace as db
import sys
from glob import glob
import time
from datetime import date
from datetime import timedelta
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



# Parameters
p0 = np.array([0.001, 0.001, 0.001, 0.001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.1, 0.00228, 0.016]) # Second to last is
bb = np.array([1e14, 1e13, 1e12, 1e11, 1e15, 1e14, 1e13, 1e12, 1e11, 1e14, 1.8e13, -0])        # PTF 12dam, nicholl17c table 3
lab = '$P={0:.0f}$~ms, ${1:s}}}$~G'
rr = 1.2e6
mm = 1.4*Msun
tt = np.logspace(-1,3,1000)
pc = 0.0336 # period Crab
pd = 4.2e-13 # pdot, buhler14
cs = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', 'k']
ci = [0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 3, 10]
alpha = [1., 0.8, 0.6, 0.4, 1., 0.8, 0.6, 0.4, 0.2, 1., 1., 1.]


# Crab numbers
ii = 2*mm*rr**2/5
tmp = np.sqrt(3*cc**3*ii/(8*np.pi**2*rr**6)*(pc*pd))
bb[-1] = tmp

# margalit18, "If the magnetic dipole and rotational axes are aligned"
# Erot = 0.5*ii*4*np.pi**2*1e6
# print(Erot)
# tau = 0.5*3*cc**3*ii/(8*np.pi**2*rr**6)*1.e0**2/1e12**2
# trot = 4.7*3600*24
# trot = 4.1e14*1e2**-2*1e-3**2
# print(trot, trot/(3600*24))
# db()
# 
# Also, see kasen10, Eq. (2). They use R=10 km

# Plots
fig = plt.figure(figsize=(8, 1.3*3.75))
for jj in range(0, p0.size):
    tp = 3/40*mm*cc**3/(np.pi**2*rr**4)*p0[jj]**2/bb[jj]**2/(3600*24*365.25)
    ll = 2**5/3*np.pi**4/cc**3*rr**6/p0[jj]**4*bb[jj]**2*1/(1+tt/tp)**2
    tmp = '{0:.1e}'.format(bb[jj]).replace('e+','\\times10^{')
    plt.loglog(tt, ll, label=lab.format(p0[jj]*1e3, tmp), color=cs[ci[jj]], alpha=alpha[jj])


# Crab numbers
sd = 2**5/3*np.pi**4/cc**3*rr**6/pc**4*bb[-1]**2
check = 4*np.pi**2*ii*pd/pc**3
xl = 3.e37
print('Crab spin-down: {0:.1e} erg s-1'.format(sd))
print('Crab 0.3-10 keV: {0:.1e} erg s-1'.format(xl))


# Cosmetics
plt.axhline(sd, color='k', ls=':')
plt.axhline(xl, color='k', ls=':')
plt.axvspan(1, 100, color='#efefef')
plt.axvspan(3, 30, color='#afafaf')
plt.xlabel('$t$ (yr)')
plt.ylabel('$L_\mathrm{bol}$ (erg s$^{-1}$)')
plt.legend(ncol=2)
plt.gca().grid(True, which='both', ls=":")
plt.savefig('nss.pdf', bbox_inches='tight', pad_inches=0.1, dpi=300)
plt.show()


db()
