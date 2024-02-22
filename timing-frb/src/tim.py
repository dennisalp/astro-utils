'''
2018-05-19, Dennis Alp, dalp@kth.se

Detection of pulsation in photon counting data. Tests different test
statistics on simulated data.
'''

from __future__ import division, print_function
import os
import pdb
import sys
from glob import glob

import numpy as np
import matplotlib.pyplot as plt

def get_Pn(*, tt, ww, k2, nn):
    walk = ww*np.exp((-2*np.pi*1j*nn)*tt) # map to complex numbers
    walk = np.sum(walk) # sum all complex numbers
    walk = np.real(walk)**2+np.imag(walk)**2 # distance from origin
    walk = walk/k2 # normalize
    return walk

def rayleigh(*, tt, ww, k2):
    walk = ww*np.exp(-2*np.pi*1j*tt)
    walk = np.sum(walk)
    walk = np.real(walk)**2+np.imag(walk)**2
    walk = walk/k2
    return walk

def h_statistic(*, tt, ww, k2):
    mmax = 20
    mm = np.arange(mmax)+1
    Z2 = np.zeros(mmax)
    HH = np.zeros(mmax)
    for nn in mm:
        Z2[nn-1] = get_Pn(tt=tt, ww=ww, k2=k2, nn=nn)

    for ii in range(0, mmax):
        HH[ii] = np.sum(Z2[:ii+1])-4*(ii+1)+4
    return np.max(HH)

def search_periods(*, pmin, pmax, tt, ww, k2, TT, sparse):
    top_min = 0. # lowest statistic in top list
    top = np.zeros((1000000, 3))
    dis = np.zeros(10000)
    jj = 0
    while pmin < pmax:
        pdot = [0.]
        pdotmax = np.min((SPIN_DOWN_LIM*pmin**3, AGE_LIM*pmin))

        while pdot[-1] < pdotmax:
            pdot.append(pdot[-1]+0.1*pmin**2/TT**2)
        pdot = pdot[:-1]

        for ii in range(-len(pdot)+1, len(pdot)):
            pdot = 0.1*ii*pmin**2/TT**2
            ttmp = tt/pmin-0.5*pdot*tt**2/pmin**2
            HH = rayleigh(tt=ttmp, ww=ww, k2=k2)
#            HH = 0.
            dis[np.min((int(HH*100), 9999))] += 1
        
#            if HH > top_min:
#                idx = np.argmax(top[:, 2] == top_min)
#                top[idx, 0] = pmin
#                top[idx, 1] = pdot
#                top[idx, 2] = HH
#                top_min = np.min(top[:, 2])
            top[jj, :] = pmin, pdot, HH
            jj += 1

            if (np.mod(jj, 1000)==0): print('{:19.14e} {:21.14e} {:21.14e}'.format(pmin, pdot, HH))

        pmin += sparse*pmin**2/(TT*10)
    return dis, top

def load_data(*, fpath):
    dat = np.loadtxt(fpath)
    tt = dat[:, 0]
    ww = dat[:, 1]
    return tt, ww

########
# Parameters
TREF = 0.
AGE_LIM = 5e-119 # 3 years
SPIN_DOWN_LIM = 6.042267201371935e-118 # 1e6 Lsun

########
# Main
fpath = sys.argv[1]
pmin = 0.03
pmax = 1.e2
sparse = 10.

tt, ww = load_data(fpath=fpath)
tt = tt-TREF
k2 = 0.5*np.sum(ww**2)
TT = tt[-1]-tt[0]

dis, top = search_periods(pmin=pmin, pmax=pmax, tt=tt, ww=ww, k2=k2, TT=TT, sparse=sparse)
np.savetxt('out/dis_py.txt', dis.astype(np.int), fmt='%i')
nn = np.argmax(top[:,0]==0.)
top = top[:nn,:]
np.savetxt('out/top_py.txt', top, fmt='%24.17e')

pdb.set_trace()
