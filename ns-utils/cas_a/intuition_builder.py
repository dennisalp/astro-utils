'''
2020-05-26, Dennis Alp, dalp@kth.se

Create some fake events and try out some formulae.
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

p0 = 0.1
tt = 100000
v0 = 287e5
a0 = 30000 #30

truth = np.arange(0, tt, p0)
n0 = 1/p0
nd = -a0*n0/cc
pd = a0*p0/cc

# Start from the same point with initial velocity v0 and acceleration a0 away from observer
obs = truth+(v0*truth+a0*truth**2/2)/cc
np.random.seed(22)
fake = np.sort(truth[np.random.choice(truth.size, 1000, replace=False)])
np.savetxt('/Users/silver/box/phd/pro/nss/cas/dat/fake0.raw', fake)
fake = np.sort(obs[np.random.choice(obs.size, 1000, replace=False)])
np.savetxt('/Users/silver/box/phd/pro/nss/cas/dat/fake1.raw', fake)

# ransom02 eq. (32)
# SHOULD BE A FACTOR OF 0.5 NOT 2
tp = obs+0.5*nd/n0*obs**2
full = -cc*(1+v0/cc)/a0+cc*np.sqrt((1+v0/cc)**2+2*a0/cc*obs)/a0 # all orders
nov0 = -cc/a0+cc*np.sqrt(1+2*a0/cc*obs)/a0
# clark17
ph = 2*np.pi*(n0*obs+0.5*nd*obs**2)

# temporal evolution of nu
nu_test = n0-a0*truth*n0/cc

# lorimer04 eq. (6.16), the observed period as function of observed time
po = p0*(1+a0*obs/cc)

db()

'''
Why are the corrections of the time stamps to second order but 
the corrections to P and nu only first order?

obs is when we see pulses.
po is the observed period, i.e. difference between two pulses. 
The difference between two pulses is constant for constant v.
This means that it basically is the derivative of the mapping used for obs.
That means that acceleration is the first-order correction for P and nu.

This is appears wierd because the events do not constitute a continuous function.
The period is not equivalent since it is the "continuous difference" between the events.
This difference, which is not a differential, makes it of 
"a lower order", even though the units are the same.
'''
