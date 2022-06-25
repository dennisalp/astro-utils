'''
2020-05-13, Dennis Alp, dalp@kth.se
python -W ignore exp.py

The high background is produced by protons in the Earth’s
magnetosphere and depends on the satellite altitude and level of solar
activity (Carter & Read 2007). These proton flares affect 30 to 40 %
of the total XMM-Newton observation time.

Request ASCII from https://www.swift.ac.uk/analysis/nhtot/index.php (more copy-paste friendly)
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
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units
from astropy.time import Time

from hlp import *

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)



################################################################
xmm = read_log('/Users/silver/box/phd/pro/obs/met/dat/xmm_log.txt')
print_coo(xmm['coo'], '/Users/silver/box/phd/pro/obs/met/dat/xmm_log_coo.txt')
xmm['nh'] = read_nh('/Users/silver/box/phd/pro/obs/met/dat/xmm_log_nh.txt')
xmm['nh'][13347] = 1.54e+20 # Problem with willingale13? It returns 0 cm-2. Using HI4PI instead.
print_diag(xmm, 'xmm')

alp = read_log('/Users/silver/box/phd/pro/obs/met/dat/xmm_log_alp20.txt')
print_coo(alp['coo'], '/Users/silver/box/phd/pro/obs/met/dat/xmm_log_alp20_coo.txt')
alp['nh'] = read_nh('/Users/silver/box/phd/pro/obs/met/dat/xmm_log_alp20_nh.txt')
alp['nh'][13347] = 1.54e+20 # Problem with willingale13? It returns 0 cm-2. Using HI4PI instead.
print_diag(alp, 'alp')

ros = read_ros('/Users/silver/box/phd/pro/obs/met/dat/ros_log.txt')
print_coo(ros['coo'], '/Users/silver/box/phd/pro/obs/met/dat/ros_log_coo.txt')
ros['nh'] = read_nh('/Users/silver/box/phd/pro/obs/met/dat/ros_log_nh.txt')
print_ros(ros, 'ros')

nus = read_nus('/Users/silver/box/phd/pro/obs/met/dat/nus_log.txt')
print_coo(nus['coo'], '/Users/silver/box/phd/pro/obs/met/dat/nus_log_coo.txt')
nus['nh'] = read_nh('/Users/silver/box/phd/pro/obs/met/dat/nus_log_nh.txt')
print_nus(nus, 'nus')
