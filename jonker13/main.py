'''
2024-02-23, Dennis Alp, dalp@kth.se

Search for pulsation in XRT 000519 (Jonker et al. 2013).
'''

from pdb import set_trace as st
import os

import numpy as np
from scipy.optimize import minimize
from scipy.stats import poisson
import pandas as pd
import matplotlib.pyplot as plt

from astropy.io import fits

from timing_lat import search, plt_stat



def plt_img(xx, yy, rr):
    img, _, _ = np.histogram2d(yy, xx, (int(rr), int(rr)))
    plt.figure()
    plt.imshow(np.log10(img), origin='lower')


def get_poi_logl(params, xx, dat):
    aa, bb, cc = params
    fit =  aa * (xx+cc)**bb
    return np.sum(-poisson.logpmf(dat, fit))


def simulate_data(xx, params, nn):
    aa, bb, cc = params
    fit = aa * (xx+cc)**bb
    fit = np.tile(fit, (nn,1))
    dat_sim = np.random.poisson(fit)
    ll_sim = np.sum(-poisson.logpmf(dat_sim, fit), axis=1)
    return dat_sim, ll_sim



######################################################
# Parameters
ff = 'acisf00803N005_evt2.fits'
x0 = 3209.7287
y0 = 5452.7807

# jonker13 source region
rr, i0 = 60.975609, 2100

i0_fit, i1_fit = 2545, 2575
p0 = [10, -2, 1]
n_sim = 100000


################################################################
# Data processing
if os.path.exists(ff):
    dat = fits.open(ff)[1].data

    xx = dat['X']
    yy = dat['Y']

    ii = (xx-x0)**2 + (yy-y0)**2 < rr**2

    xx = xx[ii]
    yy = yy[ii]
    tt = dat['time'][ii]
    ee = dat['energy'][ii]
    np.savetxt('tt.csv', tt)
else:
    tt = np.loadtxt('tt.csv')

dt = np.diff(np.unique(tt)).min()
t0 = tt.min()
tt -= t0



################################################################
# Time series
b1 = np.arange(tt.min(), tt.max()+dt, dt)
y1, _ = np.histogram(tt, b1)
b2 = np.arange(tt.min(), tt.max()+dt, dt*10)
y2, _ = np.histogram(tt, b2)
b3 = np.arange(tt.min(), tt.max()+dt, dt*100)
y3, _ = np.histogram(tt, b3)

assert np.all(np.diff(tt) >= 0), 'events must be chronological'

plt.plot(b1[:-1], y1, drawstyle='steps-post')
plt.plot(b2[:-1], y2, drawstyle='steps-post')
plt.plot(b3[:-1], y3, drawstyle='steps-post')
plt.axvline(tt[i0], color='gray')
plt.grid()
plt.xlabel('$t$ (s)')
plt.ylabel('Count')



################################################################
# Search
tt = tt[i0:]
tt += np.random.uniform(0, dt, tt.size)

# for debugging and testing
# tt = np.arange(0, 10000, 12)

p_min = dt/2
p_max = 2000

pp, pn = search(tt, p_min, p_max, sparse=1)


# Plot
plt_stat(pn)

plt.figure()
plt.plot(pp, pn)
plt.xlabel('Period (s)')
plt.ylabel('Power ($P_1$)')
plt.grid()



################################################################
# Tail variability null hypothesis test
dat = y1[i0_fit:i1_fit]
xx = np.arange(dat.shape[0])

res = minimize(get_poi_logl, p0, args=(xx, dat))
params = res.x
ll_obs = res.fun

# Perform simulations
dat_sim, ll_sim = simulate_data(xx, params, n_sim)

# Calculate the probability of observing as likely or more extreme data
goodness = np.sum(ll_sim >= ll_obs) / n_sim

# Print results
print("Fitted parameters:")
print("a={0:.3f}, b={1:.3f}, c={2:.3f}".format(*params))
print("Probability of observing as likely or more extreme data: {0:.3f}".format(goodness))


# Plot goodness
plt.figure()
plt.hist(ll_sim, 100, density=True)
plt.axvline(ll_obs, color='gray')
plt.xlabel('$\log(L)$')
plt.ylabel('Density')
plt.grid()

# Plot light curve
aa, bb, cc = params
fit = aa * (xx+cc)**bb

plt.figure()
plt.plot(xx, dat, '-', label='Obs', drawstyle='steps-post')
plt.plot(xx, fit, '-', label='Fit')
plt.xlabel('$t$ (ACIS frame)')
plt.ylabel('Count')
plt.grid()
plt.legend()
plt.show()

st()
