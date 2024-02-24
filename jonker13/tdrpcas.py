'''
2024-02-23, Dennis Alp, dalp@kth.se

Time-domain resampled partially coherent acceleration search.

Acceleration search using the FFT as well as full searches also using
FFT for both circular and elliptical orbits.

The acceleration search simple searches a range of periods and periods
with derivatives.

The scan function has two modes: circular and elliptical orbits. These
scan ranges of physical parameters of orbits.

`accel` can search subsegments of a longer observation. 
x['tm'] defines the shortest segment, time minimum.
x['nm'] defines the lowest number of bins in a segment, number minimum.
'''

from __future__ import division, print_function
import os
from pdb import set_trace as st
import sys
import time

import numpy as np
from scipy.optimize import fsolve



# Constants, cgs
cc = 2.99792458e10 # cm s-1
GG = 6.67259e-8 # cm3 g-1 s-2
Msun = 1.989e33 # g
Rsun = 6.957e10 # cm



################################################################
# General help functions
def get_lc(x, tr):

    if x['1ppb']:
        tmp = np.zeros(int(np.ceil(tr.max()/x['dt'])))
        tmp[(tr//x['dt']).astype(np.int)] += 1 # assumes only 1 photon per bin
    else:
        tmp, _ = np.histogram(tr, np.arange(0, tr.max(), x['dt']))

    # ransom02
    # values equivalent to the mean of the data
    # Padding with the data mean is preferable to zero-padding since
    # zero-padding introduces low-frequency power into the Fourier
    # response.
    lc = np.ones(x['nn'])*tmp.mean()
    # discard an event or two that is accelerated out of the power-of-two array
    ii = np.min((tmp.size, lc.size))
    lc[:ii] = tmp[:ii]
    x['dc'] = lc[:ii].sum()

    return lc


def search(x, lc):
    ft = np.fft.rfft(lc)
    po = np.abs(ft)**2 / x['dc']
    nu = np.fft.rfftfreq(x['nn'], x['dt'])
    ii = np.argmax(nu > x['fm'])
    return nu[ii:], po[ii:]


def out(x, nu, po, pars):
    tmp, _ = np.histogram(po, x['dib'])
    x['dis'] += tmp

    ii = np.where(po > x['lim'])[0]
    if x['cou']+ii.size > x['buf']:
        x['out'] = np.concatenate((x['out'], np.zeros((x['buf'], 3))), axis=0)
        print('Note: extended output buffer')
    if ii.size > 0:
        i = x['cou']
        x['out'][i:i+ii.size, 0] = po[ii]
        x['out'][i:i+ii.size, 1] = nu[ii]
        x['out'][i:i+ii.size, 2:] = pars[1:]
        x['cou'] += ii.size

    if po.max() > x['top']:
        print('\nMaximum power:', po.max())
        if pars[0] == 'accel':
            print(nu[np.argmax(po)], 'Hz,', pars[1], 'cm s-2')
        elif pars[0] == 'circular':
            buf = 'Hz, Mp: {0:4.2f} Msun, Mc: {1:6.2f} Msun, Binary period: {2:8.3f} h, T0: {3:5.0f} s'
            print(nu[np.argmax(po)], buf.format(pars[1], pars[2], pars[3]/3600, pars[4]))
        elif pars[0] == 'elliptic':
            buf = 'Hz, Mp: {0:4.2f} Msun, Mc: {1:6.2f} Msun, Binary period: {2:8.3f} h, T0: {3:5.0f} s, e: {4:5.3f}, w {5:5.3f}'
            print(nu[np.argmax(po)], buf.format(pars[1], pars[2], pars[3]/3600, pars[4], pars[5], pars[6]))
        print('Single-trial probability:', np.exp(-po.max()), '\n')
        # x['top'] = po.max()


# Prints run time estimates
def timing(x, t0, ctr, frac):
    t1 = time.time()
    if ctr < x['ip'] or np.mod(ctr-x['ip']+1, x['pi']) == 0:
        tm2 = t1-t0
        tm3 = tm2*(1/frac-1)

        if tm2 > 86400:
            tm2 = 'elapsed: {0:4.1f} d'.format(tm2/86400.)
        elif tm2 > 3600:
            tm2 = 'elapsed: {0:4.1f} h'.format(tm2/3600.)
        elif tm2 > 60:
            tm2 = 'elapsed: {0:4.1f} min'.format(tm2/60.)
        else:
            tm2 = 'elapsed: {0:4.1f} s'.format(tm2)

        if tm3 > 86400:
            tm3 = ', left: {0:4.1f} d'.format(tm3/86400.)
        elif tm3 > 3600:
            tm3 = ', left: {0:4.1f} h'.format(tm3/3600.)
        elif tm3 > 60:
            tm3 = ', left: {0:4.1f} min'.format(tm3/60.)
        else:
            tm3 = ', left: {0:4.1f} s'.format(tm3)

        print('{0:6.2f}%, '.format(100*frac) + tm2 + tm3)



################################################################
# Search of accelerated sub-segments
# ng15
# eqs 3 and 4
# The constant Tsun is used to express these masses in solar units (lorimer04 eq 8.21)
def get_amax(PP, Mc):
    MM = 1.4
    return (2*np.pi/PP)**(4/3)*(GG*Msun/cc**3*Mc**3/(MM+Mc)**2)**(1/3)*cc


def resample(tt, aa):
    return tt - aa/cc * tt**2/2


def accel(x, tt):
    x['seg'] += 1
    print('\n\n\n################################################################\n', x['a1'], 'segment:', x['seg'])

    t0 = tt[0]
    tt = tt-t0
    t1 = tt[-1]
    t2 = t1 / 2

    amax = get_amax(t1*x['dt']/x['fr'], x['mc'])

    if x['ex']:
        tmax = resample(t1, -amax)
        x['nn'] = 2**int(np.ceil(np.log2(tmax/x['dt'])))
    else:
        x['nn'] = 2**int(np.ceil(np.log2(t1/x['dt'])))

    da = x['as']*x['p0']*cc/(x['nn']*x['dt'])**2*x['db']

    start = time.time()
    for ia, aa in enumerate(np.arange(-amax, amax, da)):

        tr = resample(tt, aa)
        if tr.min() < 0:
            continue

        lc = get_lc(x, tr)
        nu, po = search(x, lc)
        out(x, nu, po, ['accel', aa])
        timing(x, start, ia, (ia+1)/(2*amax/da+1))

    if t2 > x['tm']:
        ta = tt[tt< t2]
        tb = tt[tt>=t2]
        if ta.size > x['nm']: accel(x, ta)
        if tb.size > x['nm']: accel(x, tb)



################################################################
# Search of binary orbits
def fold(x, tr, pars):
    if x['1ppb']:
        tmp = np.zeros(int(np.ceil(tr[-1]/x['dt'])))
        tmp[(tr//x['dt']).astype(np.int)] += 1 # assumes only 1 photon per bin
    else:
        tmp, _ = np.histogram(tr, np.arange(0, tr[-1], x['dt']))

    lc = np.ones(x['nn'])*tmp.mean()
    # discard an event or two that is accelerated out of the power-of-two array
    ii = np.min((tmp.size, lc.size))
    lc[:ii] = tmp[:ii] 

    ft = np.fft.rfft(lc)
    po = np.abs(ft)**2 / x['dc']
    # ip = 1/np.sqrt(2)*(ft[1:-1]-ft[2:]) # Interbin power. Using 1/sqrt(2) instead of pi/4, see interbin.py, cf. ransom02
    # ip = np.abs(ip)**2/x['dc']
    nu = np.fft.rfftfreq(x['nn'], x['dt'])
    
    ii = np.argmax(nu > x['fm'])
    ft = ft[ii:]
    po = po[ii:]
    # ip = ip[ii-2:]
    nu = nu[ii:]

    out(x, nu, po, ip, pars)


def EE_hlp(ec, MM):
    EE = MM
    for ii in range(0, 5):
        EE = EE - (EE-ec*np.sin(EE)-MM)/(1-ec*np.cos(EE))
    return EE


def circular(x, tt):
    x['tot'] = x['nMp']*x['nMc']*x['nPb']*x['nTT']
    Mp = np.linspace(1.4, 2.5, x['nMp']) # mass pulsar
    Mc = np.linspace(0.1, 100, x['nMc']) # mass companion
    Pb = np.linspace(23700, 24300, x['nPb']) # orbital period, period binary
    
    for pb in Pb:
        TT = np.linspace(0, pb, x['nTT']) # epoch of periastron passage
        for T0 in TT:
            MM = 2*np.pi*(tt-T0)/pb
            for mp in Mp:
                for mc in Mc:
                    ar = (GG*Msun*(mp+mc)*pb**2/(4*np.pi**2))**(1/3.) # relative semi-major, lorimer04, 8.21
                    ap = ar*mc/(mp+mc) # pulsar semi-major, lorimer04, 8.22
                    rb = ap/cc*(np.cos(MM)+np.sin(MM))
                    tr = tt+rb
                    
                    fold(x, tr, ['circular', mp, mc, pb, T0])
                    x['ctr'] += 1
                    timing(x, x['start'], x['ctr'], x['ctr']/x['tot'])


def elliptic(x, tt):
    x['tot'] = x['nMp']*x['nMc']*x['nPb']*x['nTT']*x['nEc']*x['nWW']
    Mp = np.linspace(1.4, 2.5, x['nMp']) # mass pulsar
    Mc = np.linspace(0.1, 10, x['nMc']) # mass companion

    Pb = np.linspace(24000, 24300, x['nPb']) # orbital period, period binary
    Ec = np.linspace(1/x['nEc'], 1-1/x['nEc'], x['nEc']) # eccentricity
    WW = np.linspace(0, 2*np.pi*(1-1/x['nWW']), x['nWW']) # longitude of periastron

    for pb in Pb:
        TT = np.linspace(0, pb*(1-1/x['nTT']), x['nTT']) # epoch of periastron passage
        for T0 in TT:
            for ec in Ec:
                MM = 2*np.pi*(tt-T0)/pb # mean anomaly
                EE = EE_hlp(ec, MM) # eccentric anomaly, lorimer04, 8.20
                for mp in Mp:
                    for mc in Mc:
                        for ww in WW:
                            ar = (GG*Msun*(mp+mc)*pb**2/(4*np.pi**2))**(1/3.) # relative semi-major, lorimer04, 8.21
                            ap = ar*mc/(mp+mc) # pulsar semi-major, lorimer04, 8.22
                            rb = ap/cc*(np.cos(EE)-ec)*np.sin(ww)+ap/cc*np.sin(EE)*np.sqrt(1-ec**2)*np.cos(ww)
                            tr = tt+rb

                            fold(x, tr, ['elliptic', mp, mc, pb, T0, ec, ww])
                            x['ctr'] += 1
                            timing(x, x['start'], x['ctr'], x['ctr']/x['tot'])


def scan(x, tt):
    x['t0'] = tt[0]
    tt = tt-x['t0']
    x['nn'] = 2**int(np.ceil(np.log2(tt[-1]/x['dt'])))
    x['start'] = time.time()
    # x['dc'] = tt.size

    if x['nEc'] == 1:
        circular(x, tt)
    else:
        elliptic(x, tt)



################################################################
# For plotting results
def plt_res(dis, top, mode='period'):
    import matplotlib.pyplot as plt
    
    dis = np.loadtxt(dis)
    top = np.loadtxt(top)
    print('{0:e} trials'.format(dis.sum()))
    print('{0:e} null hypothesis probability'.format(1-(1-np.exp(-top[:,0].max()))**dis.sum()))

    xx=np.linspace(0, 100, 101)
    plt.figure()
    plt.semilogy(xx, np.exp(-xx))
    plt.semilogy(xx, np.insert(dis, 0, dis[0])/dis.sum(), drawstyle='steps')
    plt.xlabel('Power')
    plt.grid()


    if mode=='period':
        xx = 1/top[:,1]
        lab = 'Period (s)'
    else:
        xx = top[:,1]
        lab = 'Frequency (Hz)'

    plt.figure()
    plt.scatter(xx, top[:,2], c='k', s=top[:,0]-top[:,0].min()+3)
    plt.xlabel(lab)
    plt.ylabel('Acceleration (cm s$^{-2}$)')
    plt.grid()

    plt.show()
