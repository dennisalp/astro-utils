from __future__ import division, print_function
import os
import pdb
import time
import sys
from glob import glob

import numpy as np
import pandas as pd
from scipy.interpolate import griddata
from scipy.special import factorial


################################################################
# Help functions
def my_read(ii):
    print('read:', ii)
    return np.array(pd.read_csv(ff[ii], delim_whitespace=True, names=['PP', 'Pd', 'HH']))


################################################################
# Input
ff = glob('out/top_*')
#ff = ['out/top_000.txt', 'out/top_003.txt', 'out/top_006.txt']
#ff = [ff[0]]
ng = len(ff) # number of groups


################################################################
# Main
# Just load the data
start = time.time()

dat = []
for ii in range(0, ng):
    dat.append(my_read(ii))
print('Loaded:', time.time()-start)

# Find which one to use as reference
max_num = 0.
for ii in range(0, len(dat)):
    if (dat[ii].shape[0] > max_num):
        max_num = dat[ii].shape[0]
        ref = ii

# Star the stacking
res = np.zeros(dat[ref].shape)
idx = np.zeros(ng).astype(np.int)
for ii in range(0, max_num):
    if np.mod(ii, 10000) == 0: print('{:.2f}'.format(ii/max_num))
    res[ii,:] = dat[ref][ii,:]

    for jj in range(0, ng):
        if (jj == ref):
            continue

        # Find "neighboring" instead
        # Find period first
        while dat[ref][ii,0] > dat[jj][idx[jj],0]:
            idx[jj] += 1
            if idx[jj] == dat[jj].shape[0]:
                idx[jj] -= 1
                tmp = dat[jj][idx[jj],0]
                while dat[jj][idx[jj],0] == tmp:
                    idx[jj] -= 1
                idx[jj] += 1
                break

        # Find pdot
        while dat[ref][ii,1] > dat[jj][idx[jj],1]:
            idx[jj] += 1
            if idx[jj] == dat[jj].shape[0]:
                idx[jj] -= 1
                break

        if (dat[ref][ii,1] == 0.):
            while not (dat[jj][idx[jj],1] == 0.):
                idx[jj] -= 1

        res[ii, 2] += dat[jj][idx[jj],2]
        
print('Stacked:', time.time()-start)

# Save top
top = res[:,2] > 60
pdb.set_trace()
#np.save('pdc/top', res[top, :])
#np.save('pdc/res', res)

# Load
Hmax = 250
nn = 400
n2 = 20000
HH = np.linspace(0, Hmax, n2)
vv, bb = np.histogram(res[:,2], bins=np.linspace(0, Hmax, nn+1))
vv = vv/res.shape[0]*nn/Hmax

# Rebinning
rebin = np.zeros(nn*2)
hlp = np.zeros(nn*2)
rebin[::2] = vv
rebin[1::2] = vv
hlp[::2] = np.linspace(0, Hmax-Hmax/nn, nn)
hlp[1::2] = np.linspace(Hmax/nn, Hmax, nn)

# Plotting
#np.save('pdc/pdf', np.c_[hlp, rebin])
cdf = np.cumsum(vv[::-1]*Hmax/nn)[::-1]
hlp = (hlp[1::2]+hlp[:-1:2])/2
#np.save('pdc/cdf', np.c_[hlp, cdf])

print('Binned:', time.time()-start)
