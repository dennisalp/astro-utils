'''
2019-06-26, Dennis Alp, dalp@kth.se

Take spectra (from the14) and make XSPEC models.
'''

import os
from pdb import set_trace as db
import sys
from glob import glob
from astropy.io import fits

import numpy as np
import xspec
import matplotlib as plt
from heasp import *



################################################################
os.chdir('/Users/silver/box/phd/pro/14j/nus/')
mod = sys.argv[1]
dat = np.loadtxt('dat/' + mod + '.dat')
ene = dat[:,0]
#ene = ((ene[1:]+ene[:-1])/2)
with open('dat/' + mod + '.dat') as ff:
    header = ff.readline()
time_grid = np.array([tt[:-1] for tt in header.split()[2:]]).astype('double')



################################################################
# Compare alp19 with the14
if False:
    my = np.loadtxt('dat/w7_spec_20d.txt')
    import matplotlib.pyplot as plt
    plt.loglog(my[:,0]/1e3, (my[:,2]+my[:,3]+my[:,4])*(0.0512/3.5)**2*0.6/0.07)
    plt.loglog(ene, dat[:,3])
    plt.show()
    db()
    sys.exit()


################################################################
# This fixes a non-homogeneity; dd202c starts from 12 days whereas all others are from 10 days
if time_grid[0] == 12.:
    time_grid[0] = 10.

################################################################
# set table descriptors and the energy array
tab = table()
tab.setModelName(mod)
tab.setModelUnits('')
tab.setisRedshift(False)
tab.setisAdditive(True)
tab.setisError(False)
# set up the energies. note that the size is one greater
# than that of the np.array for the model fluxes
tab.setEnergies(ene)

tab.setNumIntParams(1)
tab.setNumAddParams(0)
tab_par = tableParameter()
tab_par.setName("Time")
tab_par.setInterpolationMethod(0)
tab_par.setInitialValue(20.0)
tab_par.setDelta(0.1)
tab_par.setMinimum(10.0)
tab_par.setBottom(10.0)
tab_par.setTop(100.0)
tab_par.setMaximum(100.0)
tab_par.setTabulatedValues(time_grid)

# and push it onto the vector of parameters
tab.pushParameter(tab_par)

# now set up the spectra. these are arbitrarily calculated, in a real program 
# this step would read a file or call a routine.
for ii, tt in enumerate(time_grid):
    spe = tableSpectrum()
    interpolate = 1.
    if tt == 10. and mod == 'dd202c': # This fixes a non-homogeneity; dd202c starts from 12 days whereas all others are from 10 days
        interpolate = 10./12.
        
    spe.setParameterValues(np.array([tt]))
    spe.setFlux(interpolate*dat[:-1,ii+1]*np.diff(ene))
    tab.pushSpectrum(spe)

# now write out the table.
tab.write('xsp/mod/' + mod + '.mod')
