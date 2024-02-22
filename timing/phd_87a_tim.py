from __future__ import division, print_function
import os
import pdb
import time
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.integrate import simps
from astropy.io import fits
import scipy.stats as sts

from joblib import Parallel, delayed
import multiprocessing



################################################################
# Search for phase peaks
def phase_model(tt, pp, pdot):
    return 2*np.pi*(tt/pp-0.5*pdot*tt**2/pp**2)

def plt_hlp(ii):
    plt.hist(np.mod(tt+1.5e-2, PP[ii]), 100)
    plt.show()
    
def get_PP():
#    uu = np.arange(1e-1, 1e2/3, sparse/(10*tt[-1]))
    uu = np.arange(1e-1, 1e2, sparse/(10*tt[-1]))
    return 1/uu[::-1]

def fold(ii):
    dist = np.sort(np.mod(tt, PP[ii]))/PP[ii]
    cdf = griddata(dist, cp, cp, method='linear', fill_value=1)
    cvm = simps((cdf-cp)**2, cp, even='first')
    
    if np.mod(ii, 100000) == 0:
        print(ii/PP.size)
        sys.stdout.flush()
        
    return cvm

def lat_fold(ii):
    walk = np.exp(-2*np.pi*1j/PP[ii]*tt)
    walk = np.sum(walk)
    walk = np.real(walk)**2+np.imag(walk)**2
    walk = 2*walk/tt.size

    if np.mod(ii, 100) == 1:
        elap = time.time() - start
        print('Completed:', np.round(ii/PP.size, 4), '\tHours left:', np.round(elap/3600*(1/ii*PP.size-1),4))
        sys.stdout.flush()
        
    return walk

def max_val(ii):
    walk, tmp = np.histogram(np.mod(tt, PP[ii]), 42)
    if np.mod(ii, 100) == 1:
        elap = time.time() - start
        print('Completed:', np.round(ii/PP.size, 4), '\tHours left:', np.round(elap/3600*(1/ii*PP.size-1),4))
        sys.stdout.flush()
        
    return np.max(walk)

def max_std(ii):
    walk, tmp = np.histogram(np.mod(tt, PP[ii]), 42)
    if np.mod(ii, 100) == 1:
        elap = time.time() - start
        print('Completed:', np.round(ii/PP.size, 4), '\tHours left:', np.round(elap/3600*(1/ii*PP.size-1),4))
        sys.stdout.flush()
        
    return np.std(walk)


def vis_che(stat):
    dist, bins = np.histogram(stat, 1000, normed=True)
    plt.semilogy(bins[:-1], dist)
    plt.semilogy(bins[:-1], sts.chi2.pdf(bins[:-1], 2))
    plt.show()


################################################################
# Load data
ff='/Users/silver/dat/nus/cra/10002001006/pro/nu10002001006A01_bc.evt'
ff='/Users/silver/dat/nus/cra/10013034004/pro/nu10013034004A01_bc.evt'
dat = fits.open(ff)[1].data
crab_frequency = 29.7066
times = dat['TIME']
grade = dat['GRADE']
xx = dat['X'].astype('float')
yy = dat['Y'].astype('float')
pi = dat['PI']

# Filters
reject = (grade <= 26) # Science grade
reg = [470, 530, 475, 535] # 10002001006
reg = [460, 550, 440, 520] # 10002001006
reg = [503, 514, 471, 480] # 10002001006
reg = [436, 454, 496, 512] # 10013034004
#reg = [  0, 999,   0, 999]
har = (pi > 135) # 135 is 7 keV I think
har = (pi > 210) & (pi < 1909) # 135 is 7 keV I think
har = (pi > 39) & (pi < 1909) # 135 is 7 keV I think
src = (xx >= reg[0]) & (xx < reg[1]) & (yy >= reg[2]) & (yy < reg[3])
src_cl = np.logical_and(src, reject)
src_cl = np.logical_and(src_cl, har)

## Make an image
#nx = np.max(xx)-np.min(xx)+1
#ny = np.max(yy)-np.min(yy)+1
#sx = int(-np.min(xx))
#sy = int(-np.min(yy))
#img, xedges, yedges = np.histogram2d(xx[reject], yy[reject], [nx, ny])
#plt.imshow(img[reg[0]+sx:reg[1]+sx, reg[2]+sy:reg[3]+sy], interpolation='nearest', origin='low')
#plt.show()

########
# Prep the data
tt = times[src_cl]
tt = tt - tt[0]
#tt = tt[-140000:]
sparse = 1 # For testing
PP = get_PP()
PP = np.linspace(0.997, 1.003, 100000)/crab_frequency
PP = PP[62000:64000] # 10002001006
PP = PP[57000:59000] # 10013034004
# cumulative probability
cp = np.arange(tt.size)/(tt.size-1)
xx = np.linspace(0, 1, 101)
pdb.set_trace()
########
# Main part
num_cor = 4
start = time.time()
Pn = Parallel(n_jobs=num_cor)(delayed(lat_fold)(ii) for ii in range(0, PP.size))
end = time.time()
np.save('/Users/silver/box/phd/pro/87a/nus/Pn', Pn)

print(end - start)
