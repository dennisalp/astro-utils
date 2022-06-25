'''
2019-06-20, Dennis Alp, dalp@kth.se

Fit Type Ia models to NuSTAR observations of SN 2014J.
'''

from __future__ import division, print_function
import os
from pdb import set_trace as db
import sys
from glob import glob

import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units
from threeML import *
#from threeML.io.package_data import get_path_of_data_file

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

# Constants, cgs
cc = 2.99792458e10 # cm s-1
GG = 6.67259e-8 # cm3 g-1 s-2
hh = 6.6260755e-27 # erg s
DD = 51.2 # kpc
pc = 3.086e18 # cm
kpc = 3.086e21 # cm
mpc = 3.086e24 # cm
kev2erg = 1.60218e-9 # erg keV-1
Msun = 1.989e33 # g
Lsun = 3.828e33 # erg s-1
Rsun = 6.957e10 # cm
Tsun = 5772 # K
uu = 1.660539040e-24 # g
SBc = 5.670367e-5 # erg cm-2 K-4 s-1
kB = 1.38064852e-16 # erg K-1
snj = SkyCoord('09h55m42.137s', '+69d40m25.40s', frame='fk5')
#position RA = 9h 55m 42.137(9)s, dec = +69d 40' 25.40(5)" (J2000.0) kelly14 

################################################################
def get_dat():
    src = dat_dir + '80002092006/nu80002092006A01_sr.pha'
    bkg = dat_dir + '80002092006/nu80002092006A01_bk.pha'
    rmf = dat_dir + '80002092006/nu80002092006A01_sr.rmf'
    arf = dat_dir + '80002092006/nu80002092006A01_sr.arf'
    obs06a = OGIPLike('obs06a', src, bkg, rmf, arf)
    obs06a.set_active_measurements('60.00001-79.')
    src = dat_dir + '80002092007/nu80002092007A01_sr.pha'
    bkg = dat_dir + '80002092007/nu80002092007A01_bk.pha'
    rmf = dat_dir + '80002092007/nu80002092007A01_sr.rmf'
    arf = dat_dir + '80002092007/nu80002092007A01_sr.arf'
    obs07a = OGIPLike('obs07a', src, bkg, rmf, arf)
    obs07a.set_active_measurements('60.00001-79.')
    return DataList(obs06a, obs07a)
#    tmp = [xyl, ogip]
#    return DataList(*my_plugins)
# Parse data when I want all data

def my_plt():
    dat['obs06a'].view_count_spectrum()
    plt.semilogy(np.logspace(0,4,10000),w7(np.logspace(0,4,10000)))
    plt.xlim([1,1e4])

################################################################
# Create the second instance, this time of a different type
ver = 'v01'
dat_dir = '/Users/silver/dat/nus/14j/' + ver + '/'
dat = get_dat()

spec = TemplateModel('w7')
w7 = PointSource('w7', ra=snj.ra.degree, dec=snj.dec.degree, spectral_shape=spec)
mod = Model(w7)
mod.w7.spectrum.main.w7.time = 20
mod.w7.spectrum.main.w7.time.fix = True

# model.link(model.source2.spectrum.main.composite.K_1, model.source1.spectrum.main.Powerlaw.K)
# Link two parameters with a law. The parameters of the law become free
# parameters in the fit. In this case we impose a linear relationship
# between the index of the log-parabolic spectrum and the index of the
# powerlaw in source2: index_2 = a * alpha_1 + b.
# law = Line()
# model.link(model.source2.spectrum.main.composite.index_2, model.source2.spectrum.main.composite.alpha_1, law)
# If you want to force them to be in a specific relationship,
# say index_2 = alpha_1 + 1, just fix a and b to the corresponding values,
# after the linking, like:
# model.source2.spectrum.main.composite.index_2.Line.a = 1.0
# model.source2.spectrum.main.composite.index_2.Line.a.fix = True
# model.source2.spectrum.main.composite.index_2.Line.b = 0.0
# model.source2.spectrum.main.composite.index_2.Line.b.fix = True

#jl = JointLikelihood(mod, dat)
#best_fit, likelihood = jl.fit()
#jl.results.display()
#
#fluxes = jl.results.get_point_source_flux(100 * units.keV, 1 * units.MeV)
#
#plot_point_source_spectra(jl.results, ene_min=0.1, ene_max=1e6, num_ene=500,flux_unit='erg / (cm2 s)')

mod.w7.spectrum.main.w7.K.prior = Uniform_prior(lower_bound=0, upper_bound=10)
mod.w7.spectrum.main.w7.K.prior = Log_uniform_prior(lower_bound=1e-5,upper_bound=1e1)
mod.display(complete=True)

ba = BayesianAnalysis(mod, dat)
np.seterr(divide='ignore', invalid='ignore')
db()
# This uses the emcee sampler
samples = ba.sample(n_walkers=30, burn_in=100, n_samples=1000)

ba.results.display()
fluxes_ba = ba.results.get_point_source_flux(1 * units.keV, 10 * units.MeV)
plot_point_source_spectra(ba.results, ene_min=0.1, ene_max=1e6, num_ene=500)

ba.results.corner_plot()


db()
