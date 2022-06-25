# This just extracts data from .fits files and dumps it into a .txt file


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
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units
from astropy.time import Time

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)





################################################################
# help functions
def kev2pi(kev):
    return (kev-1.6)/0.04

def pi2kev(pi):
    return 0.04*pi+1.6


def get_wcs(ff):
    dd = fits.open(ff)[1]
    wcs = WCS()
    wcs.wcs.crpix = [dd.header['TCRPX14'], dd.header['TCRPX15']]
    wcs.wcs.cdelt = [dd.header['TCDLT14'], dd.header['TCDLT15']]
    wcs.wcs.crval = [dd.header['TCRVL14'], dd.header['TCRVL15']]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return wcs


################################################################
# parameters
cow = SkyCoord('16h16m0.22418', '+22d16m04.8903s', frame='icrs') # AT 2018cow, bietenholz20, (244.00093408, 22.26802508)
files = sorted(glob('/Users/silver/dat/nus/cow/*/pro/*_bc.evt'))
wd = '/Users/silver/box/phd/pro/cow/nus/dat/'
nr = 100
ene_lim = 35 # 3 keV
ene_lim = 460 # 20 keV
ene_lim = 1500 # 61.6 keV
ene_lim = 835 # 35 keV
ene_max = 1585 # 65 keV
ene_lim = 335 # 15 keV
ene_lim = 210 # 10 keV
ene_max = 1909

ene_max = kev2pi(24)
ene_min = kev2pi(3)





################################################################
# main
cou = 0
prev = None
rr = np.arange(nr)
pp = np.empty(nr)
cp = np.empty(nr)
os.chdir(wd)

for ff in files[:]:
    print('\n\n'+ff)

    # load and label
    dat = fits.open(ff)[1].data
    wcs = get_wcs(ff)
    x0, y0 = wcs.all_world2pix(cow.ra, cow.dec, 1)
    tt = dat['TIME']
    pi = dat['PI'].astype(np.float)
    xx = dat['X'].astype(np.float)
    yy = dat['Y'].astype(np.float)

    ii = (pi > ene_min) & (pi < ene_max)
    tt = tt[ii]
    pi = pi[ii]
    xx = xx[ii]
    yy = yy[ii]
    
    for ii in rr:
        cp[ii] = ((xx-x0)**2+(yy-y0)**2 < ii**2).sum()
        pp[ii] = cp[ii]/(np.pi*ii**2)

    bkg = pp[-10:].mean()*0.8 # 0.8 to extrpolate to infinte radii
    bkg = pp[-1]
    bkg = (cp[-1]-cp[-5])/(np.pi*(rr[-1]**2-rr[-5]**2))
    sig = cp-bkg*np.pi*rr**2
    snr = sig/np.sqrt(bkg*np.pi*rr**2)
    max = snr[1:].argmax()+1
    print('This is for ordinary S/N')
    print('Radius: {3:.1f} pix ({4:.1f} arcsec), S/N: {0:.0f}, Net counts: {1:.0f}, Background: {2:.0f}'.format(snr[max], sig[max], bkg*np.pi*rr[max]**2, rr[max], rr[max]*2.458))

    print('Pulsed searches likely want contrast between pulse peak and trough, i.e. background is the entire non-pulsed fraction')
    print('Modifed S/N for 10% pulsed fraction')
    msn = 0.1*sig/np.sqrt(sig*0.9+bkg*np.pi*rr**2)
    max = msn[1:].argmax()+1
    print('Radius: {3:.1f} pix ({4:.1f} arcsec), mS/N: {0:.1f}, Pulsed counts: {1:.0f}, Non-pulsed background: {2:.0f}'.format(msn[max], 0.1*sig[max], sig[max]*0.9+bkg*np.pi*rr[max]**2, rr[max], rr[max]*2.458))

    # plt.plot(snr)
    # plt.plot(msn)
    # plt.show()
    # db()

    tt = tt[((xx-x0)**2+(yy-y0)**2 < rr[max]**2)]
    np.savetxt(ff.split('/')[-1].replace('_bc.evt', '.txt'), tt)

    if not prev is None:
        prev = np.sort(np.concatenate((tt, prev)))
        np.savetxt(ff.split('/')[-1].replace('B01_bc.evt', 'AB1.txt'), prev)
        plt.hist(tt, np.arange(tt.min(), tt.max(), 5800), histtype='step')
        plt.show()
        prev = None
    else:
        prev = tt
