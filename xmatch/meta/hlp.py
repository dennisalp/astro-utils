'''
2020-04-30, Dennis Alp, dalp@kth.se

Help functions for observation meta data.
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
# Read observations logs
def read_log(ff):
    print('\nReading:', ff)
    ff = open(ff)
    log = {}
    obs = []
    tar = []
    coo = []
    tt0 = []
    exp = []
    arc = []
    pne = []
    mod = []
    for ll in ff:
        if ll[0] == '#':
            continue

        ll = ll.split('|')
        if ll[4].strip() == '':
            continue
        obs.append(ll[1].strip())
        tar.append(ll[3].strip())
        coo.append(ll[4].strip() + ' ' + ll[5].strip())
        tt0.append(ll[6].strip())
        exp.append(float(ll[7].strip()))
        arc.append(ll[11].strip() == 'Y')
        pne.append(0 if ll[13].strip() == '' else float(ll[13]))
        mod.append(ll[14].strip() + ' FF' if 'LW' in ll[14].strip() else ll[14].strip())

    log['obs'] = np.array(obs)
    log['tar'] = np.array(tar)
    log['coo'] = SkyCoord(coo, unit=(units.hourangle, units.degree))
    log['tt0'] = Time(np.array(tt0))
    log['exp'] = np.array(exp)
    log['arc'] = np.array(arc)
    log['pne'] = np.array(pne)
    log['mod'] = np.array(mod)
    return log

def read_ros(ff):
    print('\nReading:', ff)
    ff = open(ff)
    log = {}
    obs = []
    tar = []
    coo = []
    tt0 = []
    exp = []
    for ll in ff:
        if ll[0] == '#':
            continue

        ll = ll.split('|')
        if not 'PSPC' in ll[2]:
            continue
        obs.append(ll[1].strip())
        exp.append(float(ll[3].strip()))
        coo.append(ll[4].strip() + ' ' + ll[5].strip())
        tar.append(ll[6].strip())
        tt0.append('1990-01-01 00:00:00' if ll[7].strip() == '' else ll[7].strip())

    log['obs'] = np.array(obs)
    log['tar'] = np.array(tar)
    log['coo'] = SkyCoord(coo, unit=(units.hourangle, units.degree))
    log['tt0'] = Time(np.array(tt0))
    log['exp'] = np.array(exp)
    return log

def read_nus(ff):
    print('\nReading:', ff)
    ff = open(ff)
    log = {}
    obs = []
    tar = []
    coo = []
    tt0 = []
    exp = []
    for ll in ff:
        if ll[0] == '#':
            continue

        ll = ll.split('|')
        obs.append(ll[5].strip())
        exp.append(float(ll[6].strip()))
        coo.append(ll[2].strip() + ' ' + ll[3].strip())
        tar.append(ll[1].strip())
        tt0.append(ll[4].strip())

    log['obs'] = np.array(obs)
    log['tar'] = np.array(tar)
    log['coo'] = SkyCoord(coo, unit=(units.hourangle, units.degree))
    log['tt0'] = Time(np.array(tt0))
    log['exp'] = np.array(exp)
    return log


def print_coo(coo, ff):
    print('Printing:', ff)
    ff = open(ff, 'w')
    for cc in coo:
        ff.write(cc.to_string() + '\n')
    ff.close()

def read_nh(ff):
    print('Reading:', ff)
    ff = open(ff)
    nh = []
    for ll in ff:
        if ll[0] == '#':
            continue

        ll = ll.split('|')
        nh.append(float(ll[8].strip()))
    return np.array(nh)
    

def print_diag(dat, lab):
    print('\nPrinting:', lab)
    ff = np.flatnonzero(np.core.defchararray.find(dat['mod'],'FF')!=-1)
    print('Total duration: {0:3.1f} Ms ({1:2.1f} yr)'.format(dat['exp'].sum()/1e6, dat['exp'].sum()/365.26/24/3600))
    print('Total pn exposure: {0:3.1f} Ms ({1:2.1f} yr)'.format(dat['pne'].sum()/1e6, dat['pne'].sum()/365.26/24/3600))
    print('Total pn exposure (in EFF, FF, or LW): {0:3.1f} Ms ({1:2.1f} yr)'.format(dat['pne'][ff].sum()/1e6, dat['pne'][ff].sum()/365.26/24/3600))
    print('Total pn exposure at |b| > 15: {0:3.1f} Ms ({1:2.1f} yr)'.format(dat['pne'][np.abs(dat['coo'].galactic.b) > 15*units.deg].sum()/1e6, dat['pne'][np.abs(dat['coo'].galactic.b) > 15*units.deg].sum()/365.26/24/3600))
    print('Total pn exposure at |b| > 15 (in EFF, FF, or LW): {0:3.1f} Ms ({1:2.1f} yr)'.format(dat['pne'][ff][np.abs(dat['coo'][ff].galactic.b) > 15*units.deg].sum()/1e6, dat['pne'][ff][np.abs(dat['coo'][ff].galactic.b) > 15*units.deg].sum()/365.26/24/3600))
    print('These values are not corrected for high background intervals.\nApproximately 30-40% is lost to high background (Carter & Read 2007).')
    
    bb = np.arange(0,90,1)
    be = np.empty(bb.size)
    for ii, b in enumerate(bb):
        be[ii] = dat['pne'][ff][np.abs(dat['coo'][ff].galactic.b) > b*units.deg].sum()
    
    plt.figure(figsize=(5, 3.75))
    plt.plot(bb, be/1e6)
    plt.xlabel('Latitude (b; deg)')
    plt.ylabel('pn exposure above b (Ms)')
    plt.savefig('/Users/silver/box/phd/pro/obs/met/out/' + lab + '_exp_latitude.pdf', bbox_inches='tight', pad_inches=0.1, dpi=300)

    nb = np.arange(19,23,0.1)
    nh, _ = np.histogram(np.log10(dat['nh'][ff]), bins=nb, weights=dat['pne'][ff]/1.e6)
    plt.figure(figsize=(5, 3.75))
    plt.plot(nb, np.insert(nh, 0, nh[0]), drawstyle='steps')
    plt.ylabel('pn exposure (Ms)')
    plt.xlabel('log$_{10}$($N_\mathrm{H}/$[cm$^{-2}$])')
    plt.savefig('/Users/silver/box/phd/pro/obs/met/out/' + lab + '_nh_dis.pdf', bbox_inches='tight', pad_inches=0.1, dpi=300)

    plt.figure(figsize=(5, 3.75))
    plt.plot((nb[:-1]+nb[1:])/2., nh.cumsum(), drawstyle='steps')
    plt.ylabel('Cumulative pn exposure (Ms)')
    plt.xlabel('log$_{10}$($N_\mathrm{H}/$[cm$^{-2}$])')
    plt.savefig('/Users/silver/box/phd/pro/obs/met/out/' + lab + '_nh_cum.pdf', bbox_inches='tight', pad_inches=0.1, dpi=300)

def print_ros(dat, lab):
    print('\nPrinting:', lab)
    print('Only PSPC (B and C) pointed observations.')
    print('Total exposure: {0:3.1f} Ms ({1:2.1f} yr)'.format(dat['exp'].sum()/1e6, dat['exp'].sum()/365.26/24/3600))
    print('Total exposure at |b| > 15: {0:3.1f} Ms ({1:2.1f} yr)'.format(dat['exp'][np.abs(dat['coo'].galactic.b) > 15*units.deg].sum()/1e6, dat['exp'][np.abs(dat['coo'].galactic.b) > 15*units.deg].sum()/365.26/24/3600))
    print('These values are corrected for occulted time intervals.')
    
    bb = np.arange(0,90,1)
    be = np.empty(bb.size)
    for ii, b in enumerate(bb):
        be[ii] = dat['exp'][np.abs(dat['coo'].galactic.b) > b*units.deg].sum()
    
    plt.figure(figsize=(5, 3.75))
    plt.plot(bb, be/1e6)
    plt.xlabel('Latitude (b; deg)')
    plt.ylabel('PSPC exposure above b (Ms)')
    plt.savefig('/Users/silver/box/phd/pro/obs/met/out/' + lab + '_exp_latitude.pdf', bbox_inches='tight', pad_inches=0.1, dpi=300)

    nb = np.arange(19,23,0.1)
    nh, _ = np.histogram(np.log10(dat['nh']), bins=nb, weights=dat['exp']/1.e6)
    plt.figure(figsize=(5, 3.75))
    plt.plot(nb, np.insert(nh, 0, nh[0]), drawstyle='steps')
    plt.ylabel('PSPC exposure (Ms)')
    plt.xlabel('log$_{10}$($N_\mathrm{H}/$[cm$^{-2}$])')
    plt.savefig('/Users/silver/box/phd/pro/obs/met/out/' + lab + '_nh_dis.pdf', bbox_inches='tight', pad_inches=0.1, dpi=300)

    plt.figure(figsize=(5, 3.75))
    plt.plot((nb[:-1]+nb[1:])/2., nh.cumsum(), drawstyle='steps')
    plt.ylabel('Cumulative PSPC exposure (Ms)')
    plt.xlabel('log$_{10}$($N_\mathrm{H}/$[cm$^{-2}$])')
    plt.savefig('/Users/silver/box/phd/pro/obs/met/out/' + lab + '_nh_cum.pdf', bbox_inches='tight', pad_inches=0.1, dpi=300)    

def print_nus(dat, lab):
    print('\nPrinting:', lab)
    print('Only FPMA.')
    print('Total exposure: {0:3.1f} Ms ({1:2.1f} yr)'.format(dat['exp'].sum()/1e6, dat['exp'].sum()/365.26/24/3600))
    print('Total exposure at |b| > 15: {0:3.1f} Ms ({1:2.1f} yr)'.format(dat['exp'][np.abs(dat['coo'].galactic.b) > 15*units.deg].sum()/1e6, dat['exp'][np.abs(dat['coo'].galactic.b) > 15*units.deg].sum()/365.26/24/3600))
    print('These values are corrected for occulted time intervals.')
    
    bb = np.arange(0,90,1)
    be = np.empty(bb.size)
    for ii, b in enumerate(bb):
        be[ii] = dat['exp'][np.abs(dat['coo'].galactic.b) > b*units.deg].sum()
    
    plt.figure(figsize=(5, 3.75))
    plt.plot(bb, be/1e6)
    plt.xlabel('Latitude (b; deg)')
    plt.ylabel('FPMA exposure above b (Ms)')
    plt.savefig('/Users/silver/box/phd/pro/obs/met/out/' + lab + '_exp_latitude.pdf', bbox_inches='tight', pad_inches=0.1, dpi=300)

    nb = np.arange(19,23,0.1)
    nh, _ = np.histogram(np.log10(dat['nh']), bins=nb, weights=dat['exp']/1.e6)
    plt.figure(figsize=(5, 3.75))
    plt.plot(nb, np.insert(nh, 0, nh[0]), drawstyle='steps')
    plt.ylabel('FPMA exposure (Ms)')
    plt.xlabel('log$_{10}$($N_\mathrm{H}/$[cm$^{-2}$])')
    plt.savefig('/Users/silver/box/phd/pro/obs/met/out/' + lab + '_nh_dis.pdf', bbox_inches='tight', pad_inches=0.1, dpi=300)

    plt.figure(figsize=(5, 3.75))
    plt.plot((nb[:-1]+nb[1:])/2., nh.cumsum(), drawstyle='steps')
    plt.ylabel('Cumulative FPMA exposure (Ms)')
    plt.xlabel('log$_{10}$($N_\mathrm{H}/$[cm$^{-2}$])')
    plt.savefig('/Users/silver/box/phd/pro/obs/met/out/' + lab + '_nh_cum.pdf', bbox_inches='tight', pad_inches=0.1, dpi=300)    
