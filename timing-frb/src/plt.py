from __future__ import division, print_function
import os
import pdb
import time
import sys
from glob import glob

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import griddata
from scipy.special import factorial
from scipy.stats import chi2

from joblib import Parallel, delayed
import multiprocessing

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)


################################################################
# Input
ng = 9
pdf = np.load('pdc/pdf.npy')
cdf = np.load('pdc/cdf.npy')

################################################################
# Main
# Load
Hmax = 200
n2 = 20000
HH = np.linspace(0, Hmax, n2)

# Plotting
fig = plt.figure(figsize=(5, 3.75))
plt.semilogy(pdf[:,0], pdf[:,1]+1e-99, lw=2)
plt.semilogy(HH+1/(2*n2), chi2.pdf(HH, 2*ng), ls='--', zorder=4)

plt.xlim([0, 300])
plt.ylim([3e-7, 1e-1])
# SIMULATED KK # KK = np.load('pdc/fake_kk_poi.npy')
# SIMULATED KK # plt.gca().axvline(np.percentile(KK, 100-100*0.997300203936740), c='#2ca02c', ls='--')
# SIMULATED KK # plt.gca().axvline(np.percentile(KK, 50), c='#2ca02c')
plt.ylabel('Probability')
plt.xlabel('$H_T$')
fig.savefig('/Users/silver/box/phd/pro/nss/cas/art/fig/pdf_cdf.pdf', bbox_inches='tight', pad_inches=0.1, dpi=300)
plt.show()

pdb.set_trace()
