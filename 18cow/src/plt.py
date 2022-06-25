'''
2020-06-02, Dennis Alp, dalp@kth.se

Help functions for visualizing the results.
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
def check_po(po):
    xx=np.arange(100)
    plt.plot(xx, np.exp(-xx))
    plt.hist(po, 100, normed=True, histtype='step')
    plt.gca().set_yscale('log')
    plt.show()
    
def check_dis(dis):
    xx=np.linspace(0, 100, 101)
    plt.figure()
    plt.semilogy(xx, np.exp(-xx))
    plt.semilogy(xx, np.insert(dis, 0, dis[0])/dis.sum(), drawstyle='steps')

cwd = '/Users/silver/box/phd/pro/cow/nus/out/*'
top = np.loadtxt(glob(cwd + 'top.txt')[0])
print(top.shape)
dis = np.loadtxt(glob(cwd + 'dis.txt')[0])
print('{0:e} trials'.format(dis.sum()))
print('{0:e} null hypothesis probability'.format(1-(1-np.exp(-top[:,0].max()))**dis.sum()))

plt.figure()
plt.scatter(top[:,1], np.abs(top[:,2])+1, c='k', s=top[:,0])
plt.gca().set_yscale('log')

check_dis(dis)
plt.show()

# plt.semilogy(aa/aa.sum(), drawstyle='steps')
# plt.semilogy(bb/bb.sum(), drawstyle='steps')
# plt.semilogy(cc/cc.sum(), drawstyle='steps')
# plt.semilogy(dd/dd.sum(), drawstyle='steps')
# plt.semilogy(ee/ee.sum(), drawstyle='steps')
# 1-(1-np.exp(-top[:,0].max()))**ee.sum()=0.0041603064539741386
# plt.show()

