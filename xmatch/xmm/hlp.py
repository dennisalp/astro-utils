'''
2020-04-30, Dennis Alp, dalp@kth.se

Help functions for cross-matching.
'''

from __future__ import division, print_function
import os
from pdb import set_trace as db
import sys
from glob import glob
import time
from datetime import date

import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.coordinates import search_around_sky
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units
from astropy.time import Time




################################################################
# Read observations logs or source catalogs, see notes in
# xmatch.py for requirements
def get_xmm_log():
    print('Loading the XMM Observation Log')
    ff = '/Users/silver/box/phd/pro/obs/xma/dat/xmm_obs_log.txt'
    with open(ff, 'r') as ff:
        content = ff.readlines()

    xmm_log = {}
    obs = []
    tar = []
    poi = []
    tt0 = []
    exp = []
    pub = []
    for ii, line in enumerate(content):
        if line[0] == '#':
            continue
        line = [ll.strip() for ll in line.strip().split('|')]
        if line[6] == '' or float(line[6]) <= 0.:
            print('hlp.py/get_xmm_log() WARNING This should be filtered upon download.')
            continue
        
        obs.append(line[1])
        tar.append(line[2])
        poi.append(line[3] + ' ' + line[4])
        tt0.append(Time(line[5]))
        exp.append(float(line[6]))
        if line[7] == '':
            pub.append('unk')
        else:
            pub.append(line[7])
    
    xmm_log['obs'] = np.array(obs)
    xmm_log['tar'] = np.array(tar)
    xmm_log['poi'] = SkyCoord(poi, unit=(units.hourangle, units.deg))
    xmm_log['tt0'] = np.array(tt0)
    xmm_log['exp'] = np.array(exp)
    xmm_log['pub'] = np.array(pub)
    return xmm_log



def get_xmm_cat():
    print('Loading 4XMM DR9')
    xmm_cat = {}
    dat = fits.open('/Users/silver/box/phd/pro/obs/xma/dat/4XMM_DR9cat_v1.0.fits')[1].data
    dat = dat[dat['MJD_START'].argsort()]
    xmm_cat['coo'] = SkyCoord(dat['RA'], dat['DEC'], unit=(units.deg, units.deg))
    xmm_cat['det'] = dat['DETID']
    xmm_cat['src'] = dat['SRCID']
    xmm_cat['num'] = dat['SRC_NUM']
    xmm_cat['obs'] = dat['OBS_ID']
    xmm_cat['tt0'] = Time(dat['MJD_START'], format='mjd')
    xmm_cat['fla'] = dat['SUM_FLAG']
    xmm_cat['cts'] = dat['EP_8_CTS']

    obs = fits.open('/Users/silver/box/phd/pro/obs/xma/dat/4xmmdr9_obslist.fits')[1].data
    ii = np.argsort(obs['OBS_ID'])
    so = obs['OBS_ID'][ii]
    si = np.searchsorted(so, xmm_cat['obs'])
    ci = np.take(ii, si, mode="clip")
    xmm_cat['tar'] = obs['TARGET'][ci]
    xmm_cat['exp'] = obs['PN_TEXP'][ci]
    xmm_cat['poi'] = SkyCoord(obs['RA'][ci], obs['DEC'][ci], unit=(units.deg, units.deg))
    return xmm_cat

def get_nus_log():
    print('Loading the NuSTAR Observation Log')
    ff = '/Users/silver/box/phd/pro/obs/xma/dat/nus_obs_log.txt'
    with open(ff, 'r') as ff:
        content = ff.readlines()

    nus_log = {}
    obs = []
    tar = []
    poi = []
    tt0 = []
    exp = []
    pub = []
    for ii, line in enumerate(content):
        if line[0] == '#':
            continue
        line = [ll.strip() for ll in line.strip().split('|')]
        if line[6] == '' or float(line[6]) <= 0.:
            print('hlp.py/get_nus_log() WARNING This should be filtered upon download.')
            continue
        
        obs.append(line[5])
        tar.append(line[1][:20])
        poi.append(line[2] + ' ' + line[3])
        tt0.append(Time(line[4]))
        exp.append(float(line[6]))
        if line[8] == '':
            pub.append('unk')
        else:
            pub.append(line[8])
    
    nus_log['obs'] = np.array(obs)
    nus_log['tar'] = np.array(tar)
    nus_log['poi'] = SkyCoord(poi, unit=(units.hourangle, units.deg))
    nus_log['tt0'] = np.array(tt0)
    nus_log['exp'] = np.array(exp)
    nus_log['pub'] = np.array(pub)
    return nus_log

def get_xrt_cat():
    print('Loading 2SXPS')
    xrt_cat = {}
    dat = fits.open('/Users/silver/box/phd/pro/obs/xma/dat/2SXPS_Detections.fits')[1].data
    obs = fits.open('/Users/silver/box/phd/pro/obs/xma/dat/2SXPS_Datasets.fits')[1].data
    
    obs = obs[obs['StartTime_UTC'].argsort()]
    ii = np.argsort(obs['ObsID'])
    si = np.searchsorted(obs['ObsID'], dat['ObsID'], sorter=ii)
    dat = dat[ii[si].argsort()]
    si = np.searchsorted(obs['ObsID'], dat['ObsID'], sorter=ii)
    ci = ii[si]

    xrt_cat['tt0'] = Time(obs['StartTime_UTC'][ci], scale='utc', format='iso')
    xrt_cat['exp'] = obs['ExposureUsed'][ci]
    xrt_cat['poi'] = SkyCoord(obs['RA'][ci], obs['Decl'][ci], unit=(units.deg, units.deg))
    xrt_cat['obs'] = np.array([str(x) for x in dat['ObsID']])
    ra = np.where(np.isnan(dat['RA_corrected']), dat['RA'], dat['RA_corrected'])
    de = np.where(np.isnan(dat['Decl_corrected']), dat['Decl'], dat['Decl_corrected'])
    xrt_cat['coo'] = SkyCoord(ra, de, unit=(units.deg, units.deg))
    # xrt_cat['coo'] = SkyCoord(dat['RA'], dat['Decl'], unit=(units.deg, units.deg))

    return xrt_cat


def get_cxo_cat():
    print('Loading CSC 2.0')
    cxo_cat = {}
    ff = '/Users/silver/box/phd/pro/obs/xma/dat/csccli-results.tsv'

    with open(ff, 'r') as ff:
        content = ff.readlines()

    nam = []
    ra = []
    de = []
    for ii, line in enumerate(content):
        if line[0] == '#' or line[0] == 'n':
            continue

        line = line.split('\t')
        nam.append(line[0].strip())
        ra.append(line[1].strip().replace(' ', ':'))
        de.append(line[2].strip().replace(' ', ':'))


    cxo_cat['nam'] = np.array(nam)
    cxo_cat['coo'] = SkyCoord(ra, de, unit=(units.hourangle, units.deg))
    return cxo_cat


################################################################
# Read SNe/GRB catalogs, see notes in xmatch.py for requirements
def get_sne_cat():
    print('Loading the SN catalog')
    ff = '/Users/silver/box/phd/pro/obs/xma/dat/sne_cat.txt'
    with open(ff, 'r') as ff:
        content = ff.readlines()

    sne_cat = {}
    coo = []
    tt0 = []
    hos = []
    typ = []
    nam = []
    for ii, line in enumerate(content):
        if line[0] == '#':
            continue

        coo.append(line[:32])
        if float(line[32:36]) > 1905.:
            tt0.append(Time(float(line[49:62]), format='jd'))
        else:
            tt0.append(Time('-'.join(line[32:42].split('/')), format='iso'))
        hos.append(line[85:109].strip())
        typ.append(line[109:121].strip())
        nam.append(line[136:176].strip())

    sne_cat['coo'] = SkyCoord(coo, unit=(units.hourangle, units.deg))
    sne_cat['tt0'] = np.array(tt0)
    sne_cat['hos'] = np.array(hos)
    sne_cat['typ'] = np.array(typ)
    sne_cat['nam'] = np.array(nam)
    return sne_cat

def get_grb_cat():
    print('Loading the GRB catalog')
    ff = '/Users/silver/box/phd/pro/obs/xma/dat/grb_cat.txt'
    with open(ff, 'r') as ff:
        content = ff.readlines()

    grb_cat = {}
    nam = []
    tt0 = []
    coo = []
    err = []
    enh = []
    for ii, line in enumerate(content):
        if line[0] == '#':
            continue
        if 'NO POSITION' in line:
            continue
    
        line = line.strip().split('|')
        tmp = line[0]
        nam.append(tmp.strip())
        tt0.append(Time('20'+tmp[4:6]+'-'+tmp[6:8]+'-'+tmp[8:10], format='iso'))
        coo.append(' '.join(line[1:3]).strip())
        err.append(float(line[3]))
        enh.append(line[4].strip())

    grb_cat['nam'] = np.array(nam)
    grb_cat['tt0'] = np.array(tt0)
    grb_cat['coo'] = SkyCoord(coo, unit=(units.hourangle, units.deg))
    grb_cat['err'] = np.array(err)
    grb_cat['enh'] = np.array(enh)
    return grb_cat

def get_ext_cat():
    print('Loading the extended GRB catalog (GBM+BAT). For wider crossmatches with NuSTAR')
    ff = '/Users/silver/box/phd/pro/obs/xma/dat/bat_cat.txt'
    with open(ff, 'r') as ff:
        content = ff.readlines()

    ext_cat = {}
    nam = []
    tt0 = []
    coo = []
    for ii, line in enumerate(content):
        if line[0] == '#':
            continue
        if 'NO POSITION' in line:
            continue
    
        line = line.strip().split('|')
        if line[4].strip() == 'N/A':
            continue
        
        nam.append(line[0])
        line[3] = line[3].strip()
        if not line[3][4] == '-':
            line[3] = line[3][:4] + '-' + line[3][4:6] + '-' + line[3][6:]
        tt0.append(Time(line[3].replace('T',' '), format='iso'))
        coo.append(SkyCoord(float(line[4]), float(line[5]), unit=(units.deg, units.deg)))

    ff = '/Users/silver/box/phd/pro/obs/xma/dat/gbm_cat.txt'
    with open(ff, 'r') as ff:
        content = ff.readlines()
    for ii, line in enumerate(content):
        if line[0] == '#':
            continue
        if 'NO POSITION' in line:
            continue
    
        line = line.strip().split('|')
        nam.append(line[1])
        tt0.append(Time(line[4], format='iso'))
        coo.append(SkyCoord(line[2], line[3], unit=(units.hourangle, units.deg)))
        
    ext_cat['nam'] = np.array(nam)
    ext_cat['tt0'] = np.array(tt0)
    ext_cat['coo'] = SkyCoord(np.array(coo))
    return ext_cat



################################################################
# The object header involves some string parsing related to the
# transient.
def mk_obj_header(cat, cc):
    buf = 3*'\n' + 64*'#' + '\n'

    if 'typ' in cat.keys():
        buf += 'Date: {0:<23s} Type: {1:<12s} Name: {2:<40s}\n'.format(
            cat['tt0'][cc].iso,
            cat['typ'][cc],
            cat['nam'][cc])
    else:
        buf += 'Date: {0:<23s} Name: {1:<40s}\n'.format(
            cat['tt0'][cc].iso,
            cat['nam'][cc])

    if 'hos' in cat.keys():
        buf += 'Coordinates: {0:<29s} Host: {1:<24s} (RA, dec: {2:s}; l, b: {3:s})\n'.format(
            cat['coo'][cc].to_string('hmsdms'),
            cat['hos'][cc],
            cat['coo'][cc].to_string(),
            cat['coo'][cc].galactic.to_string())
    else:
        buf += 'Coordinates: {0:<29s} (RA, dec: {1:s}; l, b: {2:s})\n'.format(
            cat['coo'][cc].to_string('hmsdms'),
            cat['coo'][cc].to_string(),
            cat['coo'][cc].galactic.to_string())
    return buf

# The text provides the information about the observation or
# catalog entrty related to the transient.
def mk_txt(log, ll, ii, dd, dt, dp, days=True):
    if 'pub' in log.keys():
        pub = log['pub'][ll]
    else:
        pub = 'yes'

    if 'tar' in log.keys():
        tar = log['tar'][ll]
    else:
        tar = 'unk'

    if not days:
        buf = 'Epoch: {2:<12s} ({6:6d} ks) Off-axis: {8:>4.1f} deg Target: {5:<20s} Obs. ID: {1:10s} Exposure: {3:3.0f} ks Pointing: {4:<29s} Public: {7:<10s}\n'
        return buf.format(
            dd,
            log['obs'][ll],
            log['tt0'][ll].iso,
            log['exp'][ll]/1e3,
            log['poi'][ll].to_string('hmsdms'),
            tar,
            dt,
            pub,
            dp[ii])
    elif dd is None:
        buf = 'Epoch: {2:<12s} ({6:6d} days) Off-axis: {8:>4.1f} arcmin Target: {5:<20s} Obs. ID: {1:10s} Exposure: {3:3.0f} ks Pointing: {4:<29s} Public: {7:<10s}\n'
    elif dt is None:
        return 'Object: {0:<22s} Offset: {1:4.1f} arcsec\n'.format(log['nam'][ll], dd[ii])
    else:
        dd = '{0:4.1f}'.format(dd[ii])
        buf = 'Epoch: {2:<12s} ({6:6d} days) Offset: {0:>6s} arcsec Off-axis: {8:>4.1f} arcmin Target: {5:<20s} Obs. ID: {1:s} Exposure: {3:3.0f} ks Pointing: {4:<29s} Public: {7:<10s}\n'

    return buf.format(
        dd,
        log['obs'][ll],
        log['tt0'][ll].iso,
        log['exp'][ll]/1e3,
        log['poi'][ll].to_string('hmsdms'),
        tar,
        dt,
        pub,
        dp[ii])


################################################################
# This takes a telescope or catalog (log) and matches with a
# database of transients (SNe/GRBs; cat).
def mk_xmatch(log, cat, lab, txt, lim):
    print('Crossmatching:', lab)
    if lim > 3*units.deg:
        print('Large off-sets, assuming temporal matching')
        ic, il, dp, _ = search_around_sky(cat['coo'], log['poi'], lim)
        dd = None
        dp = dp.deg
        xmatch_times(log, cat, lab, txt, ic, il, dd, dp)
    elif 'nam' in log.keys():
        ic, il, dd, _ = search_around_sky(cat['coo'], log['coo'], lim)
        dd = dd.arcsec
        xmatch_sources(log, cat, lab, txt, ic, il, dd)
        
    elif 'coo' in log.keys():
        ic, il, dd, _ = search_around_sky(cat['coo'], log['coo'], lim)
        dd = dd.arcsec
        dp = log['poi'][il].separation(cat['coo'][ic]).arcmin
        xmatch_detections(log, cat, lab, txt, ic, il, dd, dp)
        
    else:
        ic, il, dp, _ = search_around_sky(cat['coo'], log['poi'], lim)
        dd = None
        dp = dp.arcmin
        xmatch_detections(log, cat, lab, txt, ic, il, dd, dp)

def xmatch_times(log, cat, lab, txt, ic, il, dd, dp):
    ff = open('/Users/silver/box/phd/pro/obs/xma/out/' + lab + '.txt', 'w')
    oo = txt + '\n\Simultaneous off-axis NuSTAR 0-bounce GRB observations.'

    last = -1
    for ii in range(ic.size):
        cc = ic[ii]
        ll = il[ii]
        dt = int(np.round((log['tt0'][ll]-cat['tt0'][cc]).value*24*3.6))
        if dt > -log['exp'][ll]*2/1000 and dt < 0:
            buf = mk_txt(log, ll, ii, dd, dt, dp, False)
            if not cc == last:
                oo += mk_obj_header(cat, cc)
            oo += buf
            last = cc
    
    ff.write(oo)
    ff.close()

def xmatch_detections(log, cat, lab, txt, ic, il, dd, dp):
    ff = open('/Users/silver/box/phd/pro/obs/xma/out/' + lab + '.txt', 'w')
    fa = open('/Users/silver/box/phd/pro/obs/xma/out/' + lab + '_after.txt', 'w')
    fb = open('/Users/silver/box/phd/pro/obs/xma/out/' + lab + '_before.txt', 'w')
    fc = open('/Users/silver/box/phd/pro/obs/xma/out/' + lab + '_close.txt', 'w')

    oo = txt
    oa = txt + '\n\"After\" implies observations later than 300 days after the SN date.'
    ob = txt + '\n\"Before\" implies observations earlier than 100 days before the SN date.'
    oc = txt + '\n\"Close\" implies observations between 100 days before and 300 days after the SN date.'

    last, la, lb, lc = -1, -1, -1, -1
    for ii in range(ic.size):
        cc = ic[ii]
        ll = il[ii]
        dt = int(np.round((log['tt0'][ll]-cat['tt0'][cc]).value))
        buf = mk_txt(log, ll, ii, dd, dt, dp)

        if dt < -100: # before
            if not cc == lb:
                ob += mk_obj_header(cat, cc)
            ob += buf
            lb = cc

        elif dt < 300: # close
            if not cc == lc:
                oc += mk_obj_header(cat, cc)
            oc += buf
            lc = cc

        else: # after
            if not cc == la:
                oa += mk_obj_header(cat, cc)
            oa += buf
            la = cc

        if not cc == last:
            oo += mk_obj_header(cat, cc)
        oo += buf
        last = cc
        

    ff.write(oo)
    fa.write(oa)
    fb.write(ob)
    fc.write(oc)
    ff.close()
    fa.close()
    fb.close()
    fc.close()

def xmatch_sources(log, cat, lab, txt, ic, il, dd):
    ff = open('/Users/silver/box/phd/pro/obs/xma/out/' + lab + '.txt', 'w')
    oo = txt

    last = -1
    for ii in range(ic.size):
        cc = ic[ii]
        ll = il[ii]
        buf = mk_txt(log, ll, ii, dd, None, None)

        if not cc == last:
            oo += mk_obj_header(cat, cc)
        oo += buf
        last = cc
        

    ff.write(oo)
    ff.close()
