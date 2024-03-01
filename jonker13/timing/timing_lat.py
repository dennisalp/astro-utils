'''
2018-05-19, Dennis Alp, dalp@kth.se

Detection of pulsation in photon counting data. 
Follows the Fermi/LAT pulsar search Einstein@Home (Clark et al. 2017).
'''

from __future__ import division, print_function
from pdb import set_trace as st
import time
import sys

import numpy as np
from scipy.interpolate import griddata
from scipy.integrate import simps
import scipy.stats as sts
import matplotlib.pyplot as plt

from joblib import Parallel, delayed



################################################################


# Search for phase peaks
# def phase_model(tt, pp, pdot):
#     return 2*np.pi*(tt/pp-0.5*pdot*tt**2/pp**2)


# def help_cvm():
#     cumulative probability (cramer-von mises)
#     cp = np.arange(tt.size)/(tt.size-1)
#     xx = np.linspace(0, 1, 101)
#     return cp, xx


# def fold_cvm(ii):
#     dist = np.sort(np.mod(tt, pp[ii]))/pp[ii]
#     cdf = griddata(dist, cp, cp, method='linear', fill_value=1)
#     cvm = simps((cdf-cp)**2, cp, even='first')
    
#     if np.mod(ii, 100000) == 0:
#         print(ii/pp.size)
#         sys.stdout.flush()
        
#     return cvm


# def max_val(ii):
#     walk, tmp = np.histogram(np.mod(tt, pp[ii]), 42)
#     if np.mod(ii, 100) == 1:
#         elap = time.time() - start
#         print('Completed:', np.round(ii/pp.size, 4), '\tHours left:', np.round(elap/3600*(1/ii*pp.size-1),4))
#         sys.stdout.flush()
        
#     return np.max(walk)


# def max_std(ii):
#     walk, tmp = np.histogram(np.mod(tt, pp[ii]), 42)
#     if np.mod(ii, 100) == 1:
#         elap = time.time() - start
#         print('Completed:', np.round(ii/pp.size, 4), '\tHours left:', np.round(elap/3600*(1/ii*pp.size-1),4))
#         sys.stdout.flush()
        
#     return np.std(walk)


def get_pp(p_min, p_max, t_range, oversampling=10, sparse=1):
    uu = np.arange(1/p_max, 1/p_min, sparse/(oversampling*t_range))
    return 1/uu[::-1]


def fold_lat(tt, pp, ii, t0):
    walk = np.exp(-2*np.pi*1j/pp[ii]*tt)
    walk = np.sum(walk)
    walk = np.real(walk)**2+np.imag(walk)**2
    walk = 2*walk/tt.size

    if np.mod(ii, 10000) == 1:
        elap = time.time() - t0
        frac = ii/pp.size
        left = elap/60 * (1/ii*pp.size-1)
        print('Completed: {0:9.3f} Minutes left: {1:9.2f}'.format(frac, left))
        sys.stdout.flush()
        
    return walk


def search(tt, p_min, p_max, sparse=1):

    t_range = tt.max() - tt.min()
    pp = get_pp(p_min, p_max, t_range, sparse=sparse)
    num_cor = 4

    t0 = time.time()
    pn = Parallel(n_jobs=num_cor)(delayed(fold_lat)(tt, pp, ii, t0) for ii in range(0, pp.size))
    t1 = time.time()
    print(t1 - t0)

    return pp, np.array(pn)


def plt_stat(pn):
    dist, bins = np.histogram(pn, 1000, normed=True)
    plt.figure()
    plt.semilogy(bins[:-1], dist+dist[-1]/10)
    plt.semilogy(bins[:-1], sts.chi2.pdf(bins[:-1], 2))
    plt.axvline(sts.chi2(2).ppf(1-1/pn.size), color='gray')
    plt.grid()
    plt.xlabel('$P_1$')
    plt.xlabel('Density')


if __name__ == '__main__':

    tt = np.empty(10000)  # dummy
    p_min = 1e-1
    p_max = 1e2
    search(tt, p_min, p_max)
