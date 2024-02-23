'''
2024-02-23, Dennis Alp, dalp@kth.se

Search for pulsation in XRT 000519 (Jonker et al. 2013).

This script searches a broader parameter space using FFT.
'''

from __future__ import division, print_function
import os
from pdb import set_trace as st
import sys
from glob import glob

import numpy as np

from tdrpcas import accel, scan, plt_res



################################################################
# Parameters
x = {}
x['db'] = 1  # boost factor, 1 for science, >> 1 for debug
x['ip'] = 10  # initial prints, how many iterations to print before switching to an interval of x['pi']
x['pi'] = 500  # printing interval
x['lim'] = 14  # lowest power to save to top
x['dt'] = 1.2047  # s, sampling interval
x['dt'] = 4.8186  # s, sampling interval
# (tt[-1]-tt[0])/2**14 = 1.2046456746462377 (choose 'dt' slightly higher than this such that a little padding is required)
x['fm'] = 0.0005  # Hz, minimum frequency, to discard frequencies affected by instrumental effects
x['ex'] = True  # extend the array to accomodate events boosted beyond the power-of-two array
x['1ppb'] = False  # At most 1 photon per bin? (i.e. bins >> photons)

# Accelerated
x['mode'] = 'accel'
x['mc'] = 100000000  # Msun, companion mass
x['fr'] = .99  # maximum fraction of orbit of the segment
x['tm'] = 9e9  # s, shortest segment
x['nm'] = 300  # fewest counts in a segment
x['as'] = 1.0  # acceleration step size, units of frequency bins

# Binary scan
# x['mode'] = 'scan'
# x['ctr'] = 0
# x['nMp'] = 1  # number of pulsar masses
# x['nMc'] = 30  # number of companion masses
# x['nPb'] = 1  # number of binary periods
# x['nTT'] = 20  # number of epoch of periastron passage
# x['nEc'] = 20  # number of eccentricity steps (1 for circular)
# x['nWW'] = 20  # number of longitude of periastron

# Source files
x['a1'] = 'tt_prepared.csv'
x['ext'] = '.'+x['a1'].split('.')[-1]
files = [x['a1']]

# Preparation
x['p0'] = 2*x['dt']
x['n0'] = 1/x['p0']
x['buf'] = 1000000
x['dis'] = np.zeros(100).astype(np.int)
x['dib'] = np.linspace(0, 100, 101)



################################################################
for ii, ff in enumerate(files):
    x['top'] = x['lim']
    x['cou'] = 0
    x['seg'] = 0
    tt = np.loadtxt(ff)

    if x['mode'] == 'accel':
        x['out'] = np.zeros((x['buf'], 3))
        accel(x, tt)
    elif x['nEc'] == 1:
        x['out'] = np.zeros((x['buf'], 6))
        scan(x, tt)
    else:
        x['out'] = np.zeros((x['buf'], 8))
        scan(x, tt)

    top = '_' + x['mode'] + '_top'+x['ext']
    top = ff.replace(x['ext'], top)
    dis = '_' + x['mode'] + '_dis'+x['ext']
    dis = ff.replace(x['ext'], dis)

    np.savetxt(top, x['out'][:x['cou'],:])
    np.savetxt(dis, x['dis'])


plt_res(dis, top, mode='frrequency')
