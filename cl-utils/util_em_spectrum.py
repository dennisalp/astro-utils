'''
2018 December 13, Dennis Alp, dalp@kth.se

Some basic unit conversions for the electromagnetic spectrum.
'''

from __future__ import division, print_function
import os
from pdb import set_trace as db
import sys
from glob import glob

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

# Constants, cgs
cc = 2.99792458e10 # cm s-1
hh = 6.6260755e-27 # erg s
DD = 51.2 # kpc
pc = 3.086e18 # cm
kpc = 3.086e21 # cm
mpc = 3.086e24 # cm
kev2erg = 1.60218e-9 # erg keV-1
Msun = 1.989e33 # g
Lsun = 3.828e33 # erg s-1
uu = 1.660539040e-24 # g

def nu2lam(nu):
    return cc/nu

def nu2ev(nu):
    return hh*nu/kev2erg
    
def lam2ev(lam):
    return hh*cc/lam/kev2erg

def lam2nu(lam):
    return cc/lam

def ev2lam(ev):
    return hh*cc/ev/kev2erg

def ev2nu(ev):
    return ev/hh*kev2erg
