'''
2020-10-03, Dennis Alp, dalp@kth.se

X-ray M31 catalog of stiele11. Sky images from XMM, CXO, XRT,
PanSTARRS, DSS, and HST/ACS.
'''

from __future__ import division, print_function
import os
from pdb import set_trace as db
import sys
import time
from glob import glob
from datetime import date
import subprocess

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import scipy.stats as sts
from scipy.ndimage import gaussian_filter
from astropy.coordinates import search_around_sky
from astropy.stats import LombScargle
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.time import Time
from xmm_lc import get_cxo_lc

################################################################
def get_xmm(ff):
    dat = fits.open(ff)[1].data
    rat = dat['EP_8_RATE']
    flx = dat['EP_8_FLUX']
    err = dat['EP_8_FLUX_ERR']
    ra = dat['SC_RA']
    de = dat['SC_DEC']
    coo = SkyCoord(ra, de, unit=(u.deg, u.deg))
    t0 = dat['MJD_START']
    return rat, flx, err, coo, t0

def get_cat(ff):
    dat = fits.open(ff)[1].data
    src = dat['XMMLPt']
    rem = dat['Remarks']
    cla = dat['Class']
    lum = dat['CFlux']*f2l
    rr = dat['CRate']
    ee = dat['e_CRate']

    r1 = dat['RAh']
    r2 = dat['RAm']
    r3 = dat['RAs']
    d1 = dat['DE-']
    d2 = dat['DEd']
    d3 = dat['DEm']
    d4 = dat['DEs']

    coo = []
    for ii, cc in enumerate(r1):
        buf = '{0:d}:{1:d}:{2:.2f} {3:s}{4:d}:{5:d}:{6:.1f}'
        buf = buf.format(r1[ii], r2[ii], r3[ii], d1[ii], d2[ii], d3[ii], d4[ii])
        coo.append(buf)

    coo = SkyCoord(coo, unit=(u.hourangle, u.deg))
    
    return src, cla, rem, rr, ee, coo, lum

def get_var(ff):
    dat = fits.open(ff)[1]
    dat.header['TNULL154'] = '  ---'
    dat.header['TNULL155'] = '  ---'
    src = dat.data['XMMLPt']
    svm = dat.data['svarmax']
    fvm = dat.data['fvarmax']
    return src, svm, fvm

def mk_plt(ss, ra, de, rr, cla, rem, lum, var):

    fig = plt.figure(figsize=(10, 5))
    gs = fig.add_gridspec(2, 3)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])
    ax4 = fig.add_subplot(gs[1, 0])
    ax5 = fig.add_subplot(gs[1, 1])
    ax6 = fig.add_subplot(gs[1, 2])

    axs = [ax1, ax2, ax3, ax4, ax5, ax6]
    for kk, cat in enumerate(cats):
        plt_hips(fov[kk], axs[kk], cat, ra, de, lab[kk])
    
    
    plt.subplots_adjust(wspace=0, hspace=0.2)
    plt.subplots_adjust(left=0, right=1)

    plt.suptitle(title.format(ss, ra, de, rr, lum, cla, rem))
    if var:
        out = 'var/{0:d}.pdf'
    else:
        out = 'con/{0:d}.pdf'

    plt.savefig(out.format(ss), bbox_inches='tight', pad_inches=0.1, dpi=300)
    plt.close()
    # plt.show()








def mk_plt_lc(ss, ra, de, rr, cla, rem, lum, var, tt, r2, ll, ee, cxo):
    fig = plt.figure(figsize=(10, 5))
    gs = fig.add_gridspec(3, 6)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])
    ax4 = fig.add_subplot(gs[0, 3])
    ax5 = fig.add_subplot(gs[0, 4])
    ax6 = fig.add_subplot(gs[0, 5])
    ax7 = fig.add_subplot(gs[1:, :])

    axs = [ax1, ax2, ax3, ax4, ax5, ax6]
    for kk, cat in enumerate(cats):
        plt_hips(fov[kk], axs[kk], cat, ra, de, lab[kk])
    
    
    ax7.errorbar(tt, ll, yerr=ee, fmt='o', ms=4, label='XMM/EPIC')
    ax7.errorbar(cxo['tt'], cxo['ff'], yerr=(cxo['lo'], cxo['hi']), fmt='o', ms=4, label='CXO/ACIS')
    ax7.set_yscale('log')
    ax7.axhline(1.e39, c='k')
    ax7.axhline(3.e37, c='k')
    ax7.set_xlabel('Time (MJD)')
    ax7.set_ylabel('Luminosity (erg s-1)')
    ax7.grid(True, which='both', ls=":")
    ax7.legend(loc='best')

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.subplots_adjust(left=0.06, right=1)

    avg = np.average(r2)
    plt.suptitle(title.format(ss, ra, de, rr, lum, cla, rem))
    if var:
        out = 'var/{0:d}.pdf'
    else:
        out = 'con/{0:d}.pdf'

    plt.savefig(out.format(ss), bbox_inches='tight', pad_inches=0.1, dpi=300)
    plt.close()
    # plt.show()

# https://aladin.u-strasbg.fr/hips/list
def plt_hips(fov, ax, cat, ra, de, lab):
    url = hips.format(cat, hips_width, hips_height, fov, ra, de)
    url = fits.open(url)
    if 'XRT' in cat:
        img = url[0].data
        img = np.log10(img+img.max()/100.)
        img = img/img.max()
    else:
        img = np.zeros((hips_width, hips_height, 3))
        img[:,:,0] = url[0].data[0]
        img[:,:,1] = url[0].data[1]
        img[:,:,2] = url[0].data[2]

        norm = img.max(axis=(0,1))
        if np.min(np.abs(norm)) > 1.e-12:
            img = img/norm
    
    ax.imshow(img, origin='lower', interpolation='nearest')
    ax.set_title(lab)
    ax.plot()
    xx = hips_width/2
    xs = hips_width/10
    ax.plot([xx+0.2*xs, xx+1.2*xs], [xx, xx], color=contrast_col, lw=0.8)
    ax.plot([xx, xx], [xx+0.2*xs, xx+1.2*xs], color=contrast_col, lw=0.8)

    ax.set_axis_off()











################################################################
# Swift/XRT rate to flux. From WEBPIMMS for Gamma=2 and N_H = 3.e20 cm-2
f2l = 4.*np.pi*(0.78*3.086e24)**2
title = '{0:d}: {1:10.7f}, {2:11.7f}   {3:.5f} cts s-1 ({4:.1e} erg s-1), {5:s}, {6:s}'
contrast_col = '#ff0000'

# Just read the data
xrat, xflx, xerr, xcoo, xt0 = get_xmm('/Users/silver/box/sci/lib/w/webb20/4XMM_DR9cat_v1.0.fits')
src, cla, rem, rr, ee, coo, lum = get_cat('/Users/silver/box/sci/lib/s/stiele11/table5.fits')
sr2, svm, fvm = get_var('/Users/silver/box/sci/lib/s/stiele11/table8.fits')


cwd = '/Users/silver/box/sci/doc/proposals/20_xmm_ulx/out/stiele11/'
os.chdir(cwd)

hips = 'http://alasky.u-strasbg.fr/hips-image-services/hips2fits?hips={}&width={}&height={}&fov={}&projection=TAN&coordsys=icrs&ra={}&dec={}'
img_size = 128 # pixel
binning = 1 # pixel, how much to bin the zoomed image (80 is default, results in a pixel size of 4 arcsec)
hips_width = 256
hips_height = 256
fov = np.array([4., 0.6, 4., 0.3, 3., 0.3])/60

# https://aladin.u-strasbg.fr/hips/list
cats = ['ESAVO/P/XMM/EPIC-RGB',
        'cxc.harvard.edu/P/cda/hips/allsky/rgb',
        'nasa.heasarc/P/Swift/XRT/exp',
        'CDS/P/PanSTARRS/DR1/color-z-zg-g',
        'CDS/P/DSS2/color',
        'ESAVO/P/HST/ACS-blue']
lab = ['XMM', 'CXO', 'XRT', 'PS1', 'DSS', 'HST']


################################################################
for ii, ss in enumerate(src):

    i2 = np.where(ss==sr2)[0]
    print(ii, cla[ii], rem[ii])
    if ii < int(sys.argv[1]):
        continue
    elif rr[ii] < 0.01 and not 'SNR' in cla[ii]:
        continue

    # XMM
    xmm_i = np.where(coo[ii].separation(xcoo) < 7*u.arcsec)[0]
    
    if xmm_i.size > 0: # 4XMM-DR9
        tt = np.empty(xmm_i.size)
        r2 = np.empty(xmm_i.size)
        ll = np.empty(xmm_i.size)
        ee = np.empty(xmm_i.size)
        for jj, xi in enumerate(xmm_i):
            tt[jj] = xt0[xi]
            r2[jj] = xrat[xi]
            ll[jj] = xflx[xi]*f2l
            ee[jj] = xerr[xi]*f2l

    cxo = get_cxo_lc(coo[ii])
    cxo['ff'] *= f2l
    cxo['lo'] *= f2l
    cxo['hi'] *= f2l
    
    if xmm_i.size > 0:
        if i2.size > 0 and svm[i2[0]] > 3:
            var = True
        else:
            var = False
        mk_plt_lc(ss, coo[ii].ra.deg, coo[ii].dec.deg, rr[ii], cla[ii], rem[ii], lum[ii], var, tt, r2, ll, ee, cxo)
    else:
        if i2.size > 0 and svm[i2[0]] > 3:
            var = True
        else:
            var = False
        mk_plt(ss, coo[ii].ra.deg, coo[ii].dec.deg, rr[ii], cla[ii], rem[ii], lum[ii], var)
        
db()
