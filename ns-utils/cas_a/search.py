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

from tdrpcas import accel


################################################################
x = {}
x['db'] = 1 # boost factor, 1 for science, >> 1 for debug
x['dt'] = 5.e-4 # s, sampling interval
x['fm'] = 0.01 # Hz, minimum frequency, to discard frequencies affected by instrumental effects

x['mc'] = 10 # Msun, companion mass
x['fr'] = 0.1 # maximum fraction of orbit of the segment

x['tm'] = 1000 # s, shortest segment
x['nm'] = 100 # fewest counts in a segment
x['as'] = 8 # acceleration step size, units of frequency bins

x['ip'] = 10 # initial prints
x['pi'] = 100 # printing interval
x['lim'] = 18

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
    x['out'] = np.zeros((x['buf'], 3))
    tt = np.loadtxt(ff)

    accel(x, tt)
    np.savetxt(ff.replace('/dat/', '/out/').replace('.txt', '_top.txt'), x['out'][:x['cou'],:])
    np.savetxt(ff.replace('/dat/', '/out/').replace('.txt', '_dis.txt'), x['dis'])
