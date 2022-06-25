# This just extracts data from .fits files and dumps it into a .npy file

from __future__ import division, print_function
import os
from pdb import set_trace as db
import sys
from glob import glob

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.ndimage import gaussian_filter
from scipy.interpolate import griddata
from scipy.optimize import curve_fit
import scipy.stats as sts

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

# parameters
files = glob('/Users/silver/dat/lat/87a/evt_flt_bc_??.fits')
wd = '/Users/silver/box/phd/pro/87a/lat/dat/'
src_rad = 1
sig = 10
subnx = 32
subny = 32
subx2 = 256
suby2 = 256
ene_min = 100
ene_max = 300000
coo = SkyCoord('05h35m27.9875s', '-69d16m11.107s', frame='icrs')

# alloc
out = np.zeros((9999999, 4))
os.chdir(wd)

################################################################
# help functions
def help_gauss(dummy, x0, y0, aa, sig):
#    print(x0, y0, aa, sig)
    xx = np.arange(0,subnx)
    yy = np.arange(0,subny)[:,None]
    xx = xx-x0
    yy = yy-y0
    rr = xx**2+yy**2
    return np.reshape(aa*np.exp(-0.5*rr/sig**2), subnx*subny)

def help_gauss_const(dummy, x0, y0, aa, sig, aa2, sig2, cc):
    print(x0, y0, aa, sig, aa2, sig2, cc)
    xx = np.arange(0, nx)
    yy = np.arange(0, ny)[:,None]
    xx = xx-x0
    yy = yy-y0
    rr = xx**2+yy**2
    rr = np.where(rr > 30**2, np.inf, rr)
    global last
    last = aa*np.exp(-0.5*rr/sig**2)+aa2*np.exp(-0.5*rr/sig2**2)+cc
    return np.reshape(last, nx*ny)

def radial_profile(dat, wei, bkg):
    x0 = subx2//2
    y0 = suby2//2
    xx, yy = np.indices((dat.shape))
    rr = np.sqrt((xx - x0)**2 + (yy - y0)**2)

    nn = 1000
    steps = np.linspace(1, subx2//2, nn)
    snr = np.zeros(nn)
    src = np.zeros(nn)
    
#   the weights makes the tail of the radial profile flat
    for ii, step in enumerate(steps):
        ind = rr < step
        src[ii] = np.sum((dat[ind]-bkg)*wei[ind])
        snr[ii] = src[ii]/np.sqrt(np.sum(dat[ind]*wei[ind]**2))
    return steps, snr

def get_Pn(tt, ww, k2, nn):
    walk = ww*np.exp((-2*np.pi*1j*nn)*tt) # map to complex numbers
    walk = np.sum(walk) # sum all complex numbers
    walk = np.real(walk)**2+np.imag(walk)**2 # distance from origin
    walk = walk/k2 # normalize
    return walk

def h_statistic(tt, ww, k2):
    mmax = 20
    mm = np.arange(mmax)+1
    Z2 = np.zeros(mmax)
    HH = np.zeros(mmax)
    for nn in mm:
        Z2[nn-1] = get_Pn(tt, ww, k2, nn)

    for ii in range(0, mmax):
        HH[ii] = np.sum(Z2[:ii+1])-4*(ii+1)+4
    return np.max(HH)


################################################################
# main
cou = 0
flux = np.empty(len(files))
ftot = np.empty(len(files)) # Use this to calibrate exposure time
for nf, ff in enumerate(files[:]):
    print(ff)

    # load and label
    dat = fits.open(ff)[1].data
    tt = dat['TIME']
    ee = dat['ENERGY']
    ra = dat['RA']
    de = dat['DEC']

    cc = SkyCoord(ra*u.degree, de*u.degree, frame='fk5', unit='deg')
    ii = (coo.separation(cc).value < src_rad)
    
#    # make image
#    nx = np.max(xx)-np.min(xx)+1
#    ny = np.max(yy)-np.min(yy)+1
#    img, xedges, yedges = np.histogram2d(xx, yy, [nx, ny])

    # save and move on
    nn = ii.sum()
    out[cou:cou+nn, :] = np.vstack((tt[ii], ee[ii], ra[ii], de[ii])).T
    cou += nn
    flux[nf] = nn
    ftot[nf] = tt.size



################################################################
# Light curve
ferr = np.sqrt(flux)
ftot = ftot/np.mean(ftot)
flux = flux/ftot
ferr = ferr/ftot
#flux[0] = flux[0]*365/150
#ferr[0] = ferr[0]*365/150
#flux[-1] = flux[-1]*365/75
#ferr[-1] = ferr[-1]*365/75
plt.errorbar(np.arange(nf+1), flux, yerr=ferr)
plt.show()
db()

################################################################
# energy filter
out = out[:cou, :]
ene_cut = (out[:, 1] > ene_min) & (out[:, 1] <= ene_max)
out = out[ene_cut, :]

# make the stacked image
nx = subx2+1
ny = suby2+1
img, xedges, yedges = np.histogram2d(out[:,2], out[:,3], bins=[nx, ny], range=[[-subx2//2, subx2//2], [-suby2//2, suby2//2]])
#plt.imshow(gaussian_filter(img, 6, order=0, mode='constant', cval=0.0, truncate=4.0), origin='lower')
#plt.show()
#db()

# fit gaussian to it to get background level and weights map
guess = [nx/2, ny/2, 10, sig, 15, sig/5, 10]
pars, covar = curve_fit(help_gauss_const, np.arange(0, nx*ny), img.ravel(), guess)
bkg = pars[-1]
src = img-bkg
wei = src/bkg-1 # if high snr, empirical weights
wei = last/bkg-1 # if low snr, modeled weights

# s/n
rr, snr = radial_profile(img, wei, bkg)
#plt.plot(rr, snr)
#plt.show()
snr_max = np.argmax(snr)
snr_max = rr[snr_max]
snr_max = snr_max if snr_max < 30 else 30
print('Maximum S/R at radius:', snr_max)
plt.plot(rr, snr)

# find optimal radius for upper limit of non-detection
def dbl_gau(rr, a1, s1, a2, s2):
    fun = rr*(a1*np.exp(-0.5*rr**2/s1**2)+a2*np.exp(-0.5*rr**2/s2**2))
    return np.cumsum(fun)
rr = np.linspace(0.0001, 50, 1000)
fpsf = dbl_gau(rr, pars[2], pars[3], pars[4], pars[5])
fpsf /= fpsf[-1]
hard_bkg = 6.87918869
variance = rr**2*hard_bkg/fpsf**2
plt.loglog(rr, variance)

# spatial filter
r_lim = out[:,2]**2+out[:,3]**2 < snr_max**2
out = out[r_lim, :]

# get weights list
wei = wei[(out[:, 2]+subx2//2).astype(int), (out[:, 3]+suby2//2).astype(int)]
out = np.hstack((out, wei[:, np.newaxis]))

sub = np.sum((img-last)[nx//2-20: nx//2+20, ny//2-20: ny//2+20])
tot = np.sum((img-bkg)[nx//2-20: nx//2+20, ny//2-20: ny//2+20])
print('Residual:', np.abs(sub/tot))

# Group the observations
#cuts = np.argwhere(np.diff(out[:,0]) > 1e6)[:,0]
#cuts = np.append(cuts, out.shape[0])+1
#cuts = np.insert(cuts, 0, 0)

org = out.copy()
order = np.argsort(out[:,0])
out = out[order,:]
low_lim = np.argmin(out[:,0])
ii = 0
cuts = [0]
while (True):
    idx = np.where(out[:,0] < out[low_lim,0]+128000.)[0]
#    np.savetxt('dat/87a_g' + str(ii).zfill(2) + '.dat', out[idx, :])
    cuts.append(idx.size+cuts[-1])

    print(out[idx[0],0], out[idx[-1],0], out[idx[-1],0]-out[idx[0],0], idx.size)
    
    out[idx,0] = np.inf
    if (np.min(out[:,0]) == np.inf): break
    low_lim = np.argmin(out[:,0])
    ii += 1

cuts = np.array(cuts)
db()  
