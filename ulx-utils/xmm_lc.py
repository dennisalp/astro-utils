'''
2020-10-05, Dennis Alp, dalp@kth.se

Get an XMM light curve for a SIMBAD object. Examples:
python xmm_lc.py "SN 1987A"
python xmm_lc.py "SN 1978K"
python xmm_lc.py 123.3123123 -12.12
python xmm_lc.py 05:36:10.123 -69:35:09.107
'''

from __future__ import division, print_function
import os
from pdb import set_trace as db
import sys
import time
from glob import glob
from datetime import date
import subprocess
import re
import io

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import scipy.stats as sts
from scipy.ndimage import gaussian_filter
from astropy.coordinates import search_around_sky
from astropy.stats import LombScargle
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad
from astropy import units as u
from astropy.io import fits
from astropy.time import Time
from astropy.io import votable
from six.moves.urllib import parse, request


################################################################
def get_src():
    try:
        print('Trying to interpret input as coordinates')
        ra = sys.argv[1].replace(',', '')
        de = sys.argv[2].replace(',', '')

        if re.match(r"[0-9].[0-9]", ra):
            print('Assuming degree angles')
            return SkyCoord(ra, de, unit=(u.deg, u.deg))
        else:
            print('Assuming hour angle')
            return SkyCoord(ra, de, unit=(u.hourangle, u.deg))
    except:
        print('Failed to interpret as coordinates, searching SIMBAD')
        try:
            res = Simbad.query_object(sys.argv[1])
            ra = res['RA'][0]
            de = res['DEC'][0]
            return SkyCoord(ra, de, unit=(u.hourangle, u.deg))
        except:
            print('All interpretations failed.\nExiting')
            sys.exit()

def get_xmm(ff):
    print('Loading 4XMM-DR9')
    dat = fits.open(ff)[1].data
    rat = dat['EP_8_RATE']
    flx = dat['EP_8_FLUX']
    err = dat['EP_8_FLUX_ERR']
    ra = dat['SC_RA']
    de = dat['SC_DEC']
    coo = SkyCoord(ra, de, unit=(u.deg, u.deg))
    t0 = dat['MJD_START']
    return rat, flx, err, coo, t0

def get_xrt(ff):
    print('Loading 2SXPS detections')
    dat = fits.open(ff)[1].data
    flx = dat['Rate']*r2f # The mean count rate of this detection (i.e. in this dataset and band), corrected for pile up, vignetting, bad columns etc.
    err = (dat['Rate_neg']+dat['Rate_pos'])/2.*r2f
    ra = dat['RA']
    de = dat['Decl']
    coo = SkyCoord(ra, de, unit=(u.deg, u.deg))
    oid_xrt = dat['ObsID']
    return flx, err, coo, oid_xrt
    
def get_obs(ff):
    print('Loading 2SXPS observation log')
    dat = fits.open(ff)[1].data
    t0 = dat['StartTime_UTC']
    keep = np.where(t0!='                   ')[0]
    t0 = Time(t0[keep]).mjd
    oid_log = dat['ObsID'][keep]
    return t0, oid_log


################################################################

def get_xmm_lc(ss):
    if __name__ == "__main__":
        print('Plotting 4XMM-DR9 detections within {0:.1f} arcsec.'.format(lim))

    xmm_i = np.where(ss.separation(coo) < lim*u.arcsec)[0]
    tt = np.empty(xmm_i.size)
    rr = np.empty(xmm_i.size)
    ff = np.empty(xmm_i.size)
    ee = np.empty(xmm_i.size)
    for jj, xi in enumerate(xmm_i):
        tt[jj] = t0[xi]
        rr[jj] = rat[xi]
        ff[jj] = flx[xi]
        ee[jj] = err[xi]

    return tt, rr, ff, ee

def get_xrt_lc(ss):
    # Swift/XRT
    if __name__ == "__main__":
        print('Plotting 2SXPS detections within {0:.1f} arcsec.'.format(lim))

    xrt_i = np.where(ss.separation(scoo) < lim*u.arcsec)[0]
    idx = np.argsort(obs_src[xrt_i])
    xrt_i = xrt_i[idx]
    nn = np.unique(obs_src[xrt_i]).size
    st = np.empty(nn)
    sf = np.empty(nn)
    se = np.empty(nn)
    ww = np.zeros(nn).astype(np.int)
    
    # Stitch together the snapshots in each Obs. ID
    last = -1
    ctr = -1
    for jj, xi in enumerate(xrt_i):
        if obs_src[xi] == last:
            sf[ctr] += sflx[xi]
            se[ctr] += serr[xi]**2
        else:
            ctr += 1
            li = np.where(obs_log==obs_src[xi])[0][0]
            st[ctr] = st0[li]
            sf[ctr] = sflx[xi]
            se[ctr] = serr[xi]**2
    
        ww[ctr] += 1
        last = obs_src[xi]
    sf = sf/ww
    se = np.sqrt(se)/ww
    return st, sf, se

def get_cxo_lc(coo):
    def get_qry(coo):
        rad = lim/60.
        ra = coo.ra.deg
        de = coo.dec.deg
        dd = rad/60.
        dr = dd/np.cos(de/360.*2*np.pi)
        ra_lo = ra-dr # might lose some sources depending on how the 
        ra_hi = ra+dr # server handles wrapping RA
        de_lo = de-dd
        de_hi = de+dd
        return csc_qry.format(ra, de, ra_lo, ra_hi, de_lo, de_hi, rad)

    if __name__ == "__main__":
        print('Plotting CXO/ACIS (CSC2.0) detections within {0:.1f} arcsec.'.format(lim))
    
    qry = get_qry(coo)
    qry = {'query': qry, 'outputFormat': 'votable', 'coordFormat': 'decimal', 'nullAppearance': 'nan'}
    qry = parse.urlencode(qry).encode("utf-8")
    req = request.urlopen("http://cda.cfa.harvard.edu/csccli/getProperties", qry)
    ans = req.read()

    with io.BytesIO() as f:
        f.write(ans)
        f.seek(0)
        vot = votable.parse_single_table(f)

    cxo = {}
    cxo['ff'] = np.array(vot.array['col2'])
    cxo['lo'] = np.array(vot.array['col3'])
    cxo['hi'] = np.array(vot.array['col4'])
    cxo['lo'] = cxo['ff']-cxo['lo']
    cxo['hi'] = cxo['hi']-cxo['ff']
    cxo['tt'] = np.array(vot.array['col5'])
    return cxo

################################################################

def mk_plt(src):
    tt, rr, ff, ee = get_xmm_lc(src)
    st, sf, se = get_xrt_lc(src)
    cxo = get_cxo_lc(src)
    plt_hlp(src.ra.deg, src.dec.deg, tt, rr, ff, ee, st, sf, se, cxo)
            
def plt_hlp(ra, de, tt, rr, ff, ee, st, sf, se, cxo):
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
    
    ax7.errorbar(st, sf, yerr=se, fmt='o', ms=4, label='Swift/XRT')
    ax7.errorbar(tt, ff, yerr=ee, fmt='o', ms=4, label='XMM/EPIC')
    ax7.errorbar(cxo['tt'], cxo['ff'], yerr=(cxo['lo'], cxo['hi']), fmt='o', ms=4, label='CXO/ACIS')
    ax7.set_yscale('log')
    ax7.axhline(1.e-13, c='k')
    ax7.set_xlabel('Time (MJD)')
    ax7.set_ylabel('Flux (erg s-1 cm-2)')
    ax7.legend(loc='best')
    ax7.grid(True, which='both', ls=":")
    
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.subplots_adjust(left=0.1, right=0.96)
    
    avg = np.average(rr) if rr.size > 0 else 0.
    plt.suptitle(title.format(ra, de, avg))

# https://aladin.u-strasbg.fr/hips/list
def plt_hips(fov, ax, cat, ra, de, lab):
    url = hips.format(cat, hips_width, hips_height, fov, ra, de)
    try:
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
    except:
        img = np.ones((hips_width, hips_height))
    
    ax.imshow(img, origin='lower', interpolation='nearest')
    ax.set_title(lab)
    ax.plot()
    xx = hips_width/2
    xs = hips_width/10
    ax.plot([xx+0.2*xs, xx+1.2*xs], [xx, xx], color=contrast_col, lw=0.8)
    ax.plot([xx, xx], [xx+0.2*xs, xx+1.2*xs], color=contrast_col, lw=0.8)

    ax.set_axis_off()





################################################################
contrast_col = '#ff0000'
r2f = 4.e-11 # Swift/XRT rate to flux. From WEBPIMMS for Gamma=2 and N_H = 3.e20 cm-2
lim = 10
rat, flx, err, coo, t0 = get_xmm('/Users/silver/box/sci/lib/w/webb20/4XMM_DR9cat_v1.0.fits')
sflx, serr, scoo, obs_src = get_xrt('/Users/silver/box/sci/lib/e/evans20/2SXPS_Detections.fits')
st0, obs_log = get_obs('/Users/silver/box/sci/lib/e/evans20/2SXPS_Datasets.fits')

title = ' '.join(sys.argv[1:]) + ': {0:10.7f}, {1:11.7f}   ---  {2:.4f} cts s-1'

hips = 'http://alasky.u-strasbg.fr/hips-image-services/hips2fits?hips={}&width={}&height={}&fov={}&projection=TAN&coordsys=icrs&ra={}&dec={}'
img_size = 128 # pixel
binning = 1 # pixel, how much to bin the zoomed image (80 is default, results in a pixel size of 4 arcsec)
hips_width = 256
hips_height = 256
fov = np.array([4., 0.6, 4., 0.3, 3., 0.3])/60

# https://aladin.u-strasbg.fr/hips/list
cats = ['ESAVO/P/XMM/EPIC-RGB',
        'cxc.harvard.edu/P/cda/hips/allsky/rgb',
        'nasa.heasarc/P/Swift/XRT/exp', # 'CDS/P/SDSS9/color-alt'
        'CDS/P/PanSTARRS/DR1/color-z-zg-g',
        'CDS/P/DSS2/color',
        'ESAVO/P/HST/ACS-blue']
lab = ['XMM', 'CXO', 'XRT', 'PS1', 'DSS', 'HST']

# https://cxc.cfa.harvard.edu/csc/quick/flux.html
# https://cxc.cfa.harvard.edu/csc/cli/
# https://cxc.cfa.harvard.edu/csc/cli/adql.html
# https://cxc.cfa.harvard.edu/csc/columns/persrc.html
csc_qry = 'SELECT DISTINCT m.name,dbo.separation(m.ra,m.dec,{0:.15f},{1:.6f}) as separation,o.flux_aper_b,o.flux_aper_lolim_b,o.flux_aper_hilim_b,o.gti_mjd_obs FROM master_source m , master_stack_assoc a , observation_source o , stack_observation_assoc b , stack_source s WHERE ((( ( m.dec BETWEEN {4:.6f} AND {5:.6f} ) AND ( m.ra BETWEEN {2:.15f} AND {3:.15f} ) ) AND dbo.cone_distance(m.ra,m.dec,{0:.15f},{1:.15f})<={6:.15f})) AND (m.name = a.name) AND (s.detect_stack_id = a.detect_stack_id and s.region_id = a.region_id) AND (s.detect_stack_id = b.detect_stack_id and s.region_id = b.region_id) AND (o.obsid = b.obsid and o.obi = b.obi and o.region_id = b.region_id)ORDER BY separation ASC, name ASC'

################################################################
if __name__ == "__main__":
    src = get_src()
    mk_plt(src)
    out = '/Users/silver/Desktop/{0:.7f}_{1:.7f}.pdf'.format(src.ra.deg, src.dec.deg)
    plt.savefig(out, bbox_inches='tight', pad_inches=0.1, dpi=300)
    plt.show()
