'''
2020-10-03, Dennis Alp, dalp@kth.se

Produce light curves and get sky images of ULXs. ULX catalog of
earnshaw19b. X-ray data from 4XMM-DR9 and 2SXPS. Sky images from XMM,
CXO, XRT, PanSTARRS, DSS, and BAT.
'''

from __future__ import division, print_function
import os
from pdb import set_trace as db
import sys
import time
from glob import glob
from datetime import date
import subprocess
import io

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
from astropy.io import votable
from six.moves.urllib import parse, request

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

def get_ulx(ff):
    dat = fits.open(ff)[1].data
    src = dat['SRCID']
    rat = dat['EP_8_RATE']
    lum = dat['EP_8_LUMINOSITY']
    err = dat['EP_8_LUMINOSITY_ERR']
    dis = dat['DISTANCE_BEST'].astype(np.float64)
    ra = dat['SC_RA']
    de = dat['SC_DEC']
    coo = SkyCoord(ra, de, unit=(u.deg, u.deg))
    t0 = dat['MJD_START']
    return src, rat, lum, err, dis, coo, t0

def get_xrt(ff):
    dat = fits.open(ff)[1].data
    rat = dat['Rate']*r2f
    err = (dat['Rate_neg']+dat['Rate_pos'])/2.*r2f
    ra = dat['RA']
    de = dat['Decl']
    coo = SkyCoord(ra, de, unit=(u.deg, u.deg))
    obs_src = dat['ObsID']
    exp = dat['CorrectedExposure']
    return rat, err, coo, obs_src, exp
    
def get_obs(ff):
    dat = fits.open(ff)[1].data
    t0 = dat['StartTime_UTC']
    keep = np.where(t0!='                   ')[0]
    t0 = Time(t0[keep]).mjd
    obs_log = dat['ObsID'][keep]
    return t0, obs_log





################################################################
def get_qry(coo):
    rad = search_radius/60.
    ra = coo.ra.deg
    de = coo.dec.deg
    dd = rad/60.
    dr = dd/np.cos(de/360.*2*np.pi)
    ra_lo = ra-dr # might lose some sources depending on how the 
    ra_hi = ra+dr # server handles wrapping RA
    de_lo = de-dd
    de_hi = de+dd
    return csc_qry.format(ra, de, ra_lo, ra_hi, de_lo, de_hi, rad)







################################################################
def mk_plt_ls(ra, de, tt, rr, ll, ee, st, sr, se, nu, po):

    fig = plt.figure(figsize=(10, 5))
    gs = fig.add_gridspec(4, 6)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])
    ax4 = fig.add_subplot(gs[0, 3])
    ax5 = fig.add_subplot(gs[0, 4])
    ax6 = fig.add_subplot(gs[0, 5])
    ax7 = fig.add_subplot(gs[1:3, :])
    ax8 = fig.add_subplot(gs[3:, :])

    axs = [ax1, ax2, ax3, ax4, ax5, ax6]
    for kk, cat in enumerate(cats):
        plt_hips(fov[kk], axs[kk], cat, ra, de)
    
    
    ax7.errorbar(st, sr, yerr=se, fmt='o', ms=4, label='Swift/XRT')
    ax7.errorbar(tt, ll, yerr=ee, fmt='o', ms=4, label='XMM/EPIC')
    ax7.set_yscale('log')
    ax7.axhline(1.e39, c='k')
    ax7.axhline(3.e37, c='k')
    ax7.set_xlabel('Time (MJD)')
    ax7.set_ylabel('Luminosity (erg s-1)')
    ax7.legend(loc='best')
    ax7.grid(True, which='both', ls=":")

    ax8.plot(nu, po)
    ax8.set_xlabel('Frequency (day-1)')
    ax8.set_ylabel('L-S Power')
    
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.subplots_adjust(left=0.06, right=1)

    avg = np.average(rr)
    plt.suptitle(title.format(ra, de, avg, ii))
    if avg < 0.01:
        plt.savefig('earnshaw19b/vlow/{0:d}.pdf'.format(ii), bbox_inches='tight', pad_inches=0.1, dpi=300)
    elif avg < 0.1:
        plt.savefig('earnshaw19b/low/{0:d}.pdf'.format(ii), bbox_inches='tight', pad_inches=0.1, dpi=300)
    else:
        plt.savefig('earnshaw19b/high/{0:d}.pdf'.format(ii), bbox_inches='tight', pad_inches=0.1, dpi=300)
    plt.close()
    # plt.show()
def mk_plt(ra, de, tt, rr, ll, ee, st, sr, se, cxo):

    fig = plt.figure(figsize=(10, 5))
    gs = fig.add_gridspec(3, 6)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])
    ax4 = fig.add_subplot(gs[0, 3])
    ax5 = fig.add_subplot(gs[0, 4])
    ax6 = fig.add_subplot(gs[0, 5])
    ax7 = fig.add_subplot(gs[1:3, :])

    axs = [ax1, ax2, ax3, ax4, ax5, ax6]
    for kk, cat in enumerate(cats):
        plt_hips(fov[kk], axs[kk], cat, ra, de)
    
    
    ax7.errorbar(st, sr, yerr=se, fmt='o', ms=4, label='Swift/XRT')
    ax7.errorbar(tt, ll, yerr=ee, fmt='o', ms=4, label='XMM/EPIC')
    ax7.errorbar(cxo['tt'], cxo['ff'], yerr=(cxo['lo'], cxo['hi']), fmt='o', ms=4, label='CXO/ACIS')
    ax7.set_yscale('log')
    ax7.axhline(1.e39, c='k')
    ax7.axhline(3.e37, c='k')
    ax7.set_xlabel('Time (MJD)')
    ax7.set_ylabel('Luminosity (erg s-1)')
    ax7.legend(loc='best')
    ax7.grid(True, which='both', ls=":")

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.subplots_adjust(left=0.06, right=1)

    avg = np.average(rr)
    plt.suptitle(title.format(ra, de, avg, ii))
    if avg < 0.01:
        plt.savefig('earnshaw19b/vlow/{0:d}.pdf'.format(ii), bbox_inches='tight', pad_inches=0.1, dpi=300)
    elif avg < 0.1:
        plt.savefig('earnshaw19b/low/{0:d}.pdf'.format(ii), bbox_inches='tight', pad_inches=0.1, dpi=300)
    else:
        plt.savefig('earnshaw19b/high/{0:d}.pdf'.format(ii), bbox_inches='tight', pad_inches=0.1, dpi=300)
    plt.close()
    # plt.show()

# https://aladin.u-strasbg.fr/hips/list
def plt_hips(fov, ax, cat, ra, de):
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
    ax.set_title(cat[-13:])
    ax.plot()
    xx = hips_width/2
    xs = hips_width/10
    ax.plot([xx+0.2*xs, xx+1.2*xs], [xx, xx], color=contrast_col, lw=0.8)
    ax.plot([xx, xx], [xx+0.2*xs, xx+1.2*xs], color=contrast_col, lw=0.8)

    ax.set_axis_off()











################################################################
# Swift/XRT rate to flux. From WEBPIMMS for Gamma=2 and N_H = 3.e20 cm-2
r2f = 4.e-11
title = '{3:d}: {0:10.7f}, {1:11.7f}   ---   {2:.5f} cts s-1'
contrast_col = '#ff0000'
search_radius = 5 # arcsec, for CXO CSC

# Just read the data
xrat, xflx, xerr, xcoo, xt0 = get_xmm('/Users/silver/box/sci/lib/w/webb20/4XMM_DR9cat_v1.0.fits')
srat, serr, scoo, obs_src, exp = get_xrt('/Users/silver/box/sci/lib/e/evans20/2SXPS_Detections.fits')
st0, obs_log = get_obs('/Users/silver/box/sci/lib/e/evans20/2SXPS_Datasets.fits')
src, rat, lum, err, dis, coo, t0 = get_ulx('/Users/silver/box/sci/lib/e/earnshaw19b/catalog.fits')
cwd = '/Users/silver/box/sci/doc/proposals/20_xmm_ulx/out/'
os.chdir(cwd)

hips = 'http://alasky.u-strasbg.fr/hips-image-services/hips2fits?hips={}&width={}&height={}&fov={}&projection=TAN&coordsys=icrs&ra={}&dec={}'
img_size = 128 # pixel
binning = 1 # pixel, how much to bin the zoomed image (80 is default, results in a pixel size of 4 arcsec)
hips_width = 256
hips_height = 256
fov = np.array([4., 0.5, 4., 0.3, 3., 12.])/60

# https://aladin.u-strasbg.fr/hips/list
cats = ['ESAVO/P/XMM/EPIC-RGB',
        'cxc.harvard.edu/P/cda/hips/allsky/rgb',
        'nasa.heasarc/P/Swift/XRT/exp',
        'CDS/P/PanSTARRS/DR1/color-z-zg-g',
        'CDS/P/DSS2/color',
        'CDS/P/allWISE/color']
        # 'JAXA/P/SWIFT_BAT_FLUX']



# https://cxc.cfa.harvard.edu/csc/quick/flux.html
# https://cxc.cfa.harvard.edu/csc/cli/
# https://cxc.cfa.harvard.edu/csc/cli/adql.html
# https://cxc.cfa.harvard.edu/csc/columns/persrc.html
csc_qry = 'SELECT DISTINCT m.name,dbo.separation(m.ra,m.dec,{0:.15f},{1:.6f}) as separation,o.flux_aper_b,o.flux_aper_lolim_b,o.flux_aper_hilim_b,o.gti_mjd_obs FROM master_source m , master_stack_assoc a , observation_source o , stack_observation_assoc b , stack_source s WHERE ((( ( m.dec BETWEEN {4:.6f} AND {5:.6f} ) AND ( m.ra BETWEEN {2:.15f} AND {3:.15f} ) ) AND dbo.cone_distance(m.ra,m.dec,{0:.15f},{1:.15f})<={6:.15f})) AND (m.name = a.name) AND (s.detect_stack_id = a.detect_stack_id and s.region_id = a.region_id) AND (s.detect_stack_id = b.detect_stack_id and s.region_id = b.region_id) AND (o.obsid = b.obsid and o.obi = b.obi and o.region_id = b.region_id)ORDER BY separation ASC, name ASC'

    
################################################################
for ii, ss in enumerate(np.unique(src)):
    print(ii)
    if ii < int(sys.argv[1]):
        continue

    # XMM
    ulx_i = np.where(ss==src)[0]
    ui = ulx_i[0]
    f2l = 4*np.pi*(dis[ui].astype(np.float64)*3.086e24)**2
    xmm_i = np.where(coo[ui].separation(xcoo) < 7*u.arcsec)[0]
    
    if xmm_i.size > 0: # 4XMM-DR9
        tt = np.empty(xmm_i.size)
        rr = np.empty(xmm_i.size)
        ll = np.empty(xmm_i.size)
        ee = np.empty(xmm_i.size)
        for jj, xi in enumerate(xmm_i):
            tt[jj] = xt0[xi]
            rr[jj] = xrat[xi]
            ll[jj] = xflx[xi]*f2l
            ee[jj] = xerr[xi]*f2l
    else: # earnshaw19b, 3XMM-DR4
        tt = np.empty(ulx_i.size)
        rr = np.empty(ulx_i.size)
        ll = np.empty(ulx_i.size)
        ee = np.empty(ulx_i.size)
        for jj, ui in enumerate(ulx_i):
            tt[jj] = t0[ui]
            rr[jj] = rat[ui]
            ll[jj] = lum[ui]
            ee[jj] = err[ui]

    # CXO
    qry = get_qry(coo[ui])
    qry = {'query': qry, 'outputFormat': 'votable', 'coordFormat': 'decimal', 'nullAppearance': 'nan'}
    qry = parse.urlencode(qry).encode("utf-8")
    req = request.urlopen("http://cda.cfa.harvard.edu/csccli/getProperties", qry)
    ans = req.read()

    with io.BytesIO() as f:
        f.write(ans)
        f.seek(0)
        vot = votable.parse_single_table(f)

    cxo = {}
    cxo['ff'] = np.array(vot.array['col2'])*f2l
    cxo['lo'] = np.array(vot.array['col3'])*f2l
    cxo['hi'] = np.array(vot.array['col4'])*f2l
    cxo['lo'] = cxo['ff']-cxo['lo']
    cxo['hi'] = cxo['hi']-cxo['ff']
    cxo['tt'] = np.array(vot.array['col5'])
        
    # ct = np.empty(xmm_i.size)
    # cr = np.empty(xmm_i.size)
    # cl = np.empty(xmm_i.size)
    # ce = np.empty(xmm_i.size)
    
    # for jj, xi in enumerate(xmm_i):
    #     tt[jj] = xt0[xi]
    #     rr[jj] = xrat[xi]
    #     ll[jj] = xflx[xi]*4*np.pi*(dis[ui].astype(np.float64)*3.086e24)**2
    #     ee[jj] = xerr[xi]*4*np.pi*(dis[ui].astype(np.float64)*3.086e24)**2

            
    # Swift
    xrt_i = np.where(coo[ui].separation(scoo) < 7*u.arcsec)[0]
    idx = np.argsort(obs_src[xrt_i])
    xrt_i = xrt_i[idx]
    nn = np.unique(obs_src[xrt_i]).size
    st = np.empty(nn)
    sr = np.empty(nn)
    se = np.empty(nn)
    ww = np.zeros(nn).astype(np.int)

    # Stitch together the snapshots in each Obs. ID
    last = -1
    ctr = -1
    for jj, xi in enumerate(xrt_i):
        if obs_src[xi] == last:
            sr[ctr] += srat[xi]
            se[ctr] += serr[xi]**2
        else:
            ctr += 1
            li = np.where(obs_log==obs_src[xi])[0][0]
            st[ctr] = st0[li]
            sr[ctr] = srat[xi]
            se[ctr] = serr[xi]**2

        ww[ctr] += 1
        last = obs_src[xi]
    sr = sr/ww
    se = np.sqrt(se)/ww

    sr = sr*4*np.pi*(dis[ui].astype(np.float64)*3.086e24)**2
    se = se*4*np.pi*(dis[ui].astype(np.float64)*3.086e24)**2
    
    # ls = LombScargle(st, sr)
    # nu, po = ls.autopower()
    # plt.plot(nu, po)
    # ls.false_alarm_probability(po.max())
    mk_plt(coo[ui].ra.deg, coo[ui].dec.deg, tt, rr, ll, ee, st, sr, se, cxo)

db()
