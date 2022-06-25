'''
2019-07-19, Dennis Alp, dalp@kth.se

Try some very basic statistics to understand what happens in grppha and XSPEC.
'''

from __future__ import division, print_function
import os
from pdb import set_trace as db
import sys
from glob import glob
from datetime import date

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sts

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

cts = int(sys.argv[2])
ng = int(sys.argv[3])
goodness = np.zeros(ng)

for jj in range(ng):
    nn = float(sys.argv[1])
    nn = np.random.poisson(nn)
    evt = np.random.uniform(0, 1, nn) # The "events", the analogy would be the energy of a list of photons
    evt = np.sort(evt)
    nb = nn//cts
    edges = np.zeros(2*nb)
    pha = np.zeros(2*nb)
    for ii in range(nb):
        if ii == nb-1:
            edges[2*ii+1] = evt[(ii+1)*cts-1]
            flx = cts/(1.-edges[2*ii])
        else:
            edges[2*ii+1] = edges[2*ii+2] = evt[(ii+1)*cts-1]
            flx = cts/(edges[2*ii+1]-edges[2*ii])
        pha[2*ii] = pha[2*ii+1] = flx
    
    # Stats
    chi2 = np.sum((pha[::2]-nn)**2/(np.sqrt(cts)/cts*pha[::2])**2)
#    print('The chi2 is {0:7.2f} for {1:4d} dof, resulting in a reduced chi2 of {2:4.2f}'.format(chi2, nb, chi2/nb))
    goodness[jj] = chi2

print('Goodness:', np.sum(goodness<30.38)/ng)
print('Should be around', sts.chi2.cdf(30.38,float(sys.argv[1])/cts))
plt.hist(goodness, 100, histtype='step')
plt.axvline(30.38)
plt.show()
db()
# Plot
#plt.plot(edges, pha, 'k')
#plt.plot(edges, (1.+np.sqrt(cts)/cts)*pha, color='gray')
#plt.plot(edges, (1.-np.sqrt(cts)/cts)*pha, color='gray')
#plt.axhline(nn)
#plt.show()
#db()
