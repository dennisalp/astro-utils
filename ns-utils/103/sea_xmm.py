'''
2020-05-27, Dennis Alp, dalp@kth.se

Search for coherent pulsations.

time python -m cProfile -s cumtime search.py  
'''

from __future__ import division, print_function
import os
from pdb import set_trace as db
import sys
from glob import glob

import numpy as np

from tdrpcas import accel, scan


################################################################
x = {}
x['db'] = 1 # boost factor, 1 for science, >> 1 for debug
x['ip'] = 10 # initial prints
x['pi'] = 100 # printing interval
x['lim'] = 6
x['dt'] = 0.308881 # s, sampling interval (choose such that little padding is required)
# (tt[-1]-tt[0])/2**18=0.30888080961631204
x['fm'] = 0.002 # Hz, minimum frequency, to discard frequencies affected by instrumental effects

# Accelerated
x['mode'] = 'accel'
x['1ppb'] = False # At most 1 photon per bin? (i.e. bins >> photons)
x['mc'] = 10 # Msun, companion mass
x['fr'] = 0.1 # maximum fraction of orbit of the segment
x['tm'] = 1000 # s, shortest segment
x['nm'] = 100 # fewest counts in a segment
x['as'] = 8 # acceleration step size, units of frequency bins

# Binary scan
# x['mode'] = 'scan'
# x['1ppb'] = False # At most 1 photon per bin? (i.e. bins >> photons)
# x['ctr'] = 0
# x['nMp'] = 1 # number of pulsar masses
# x['nMc'] = 30 # number of companion masses
# x['nPb'] = 1 # number of binary periods
# x['nTT'] = 20 # number of epoch of periastron passage
# x['nEc'] = 20 # number of eccentricity steps (1 for circular)
# x['nWW'] = 20 # number of longitude of periastron

# Preparation
x['p0'] = 2*x['dt']
x['n0'] = 1/x['p0']
x['buf'] = 1000000
x['dis'] = np.zeros(100).astype(np.int)
x['dib'] = np.linspace(0, 100, 101)
x['a1'] = sys.argv[1]
files = sorted(glob('../dat/' + x['a1'] + '.txt'))



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

    np.savetxt(ff.replace('/dat/', '/out/').replace('.txt', '_' + x['mode'] + '_top.txt'), x['out'][:x['cou'],:])
    np.savetxt(ff.replace('/dat/', '/out/').replace('.txt', '_' + x['mode'] + '_dis.txt'), x['dis'])
