'''
2019-06-20, Dennis Alp, dalp@kth.se

Take spectra (from the14) and make 3ML models (astromodels).
'''

from __future__ import division, print_function
import os
from pdb import set_trace as db
import sys
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as uu

from astromodels import Band
import numpy as np
from astromodels import TemplateModelFactory
from astromodels import TemplateModel

# https://threeml.readthedocs.io/en/latest/notebooks/spectral_models.html
# 
# Template (Table) Models
# 
# 3ML (via astromodels) provides the ability to construct models in
# tabluated form. This is very useful for models that are from
# numerical simualtions. While in other software special care must be
# taken to construct table models into FITS files, in 3ML the
# construction of the table model is taken care of for you. Here is an
# example of how to build a template model from a pre-existing
# function.

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
os.chdir('/Users/silver/box/phd/pro/14j/nus/')


################################################################
mod = sys.argv[1]
dat = np.loadtxt('spectra/' + mod + '.dat')
ene = dat[:,0]
ene = ((ene[1:]+ene[:-1])/2)
with open('spectra/' + mod + '.dat') as ff:
    header = ff.readline()
time_grid = np.array([tt[:-1] for tt in header.split()[2:]]).astype('double')

# Now we define a template model factory. This takes a name, a
# description, the energy grid and an array of parameter names as
# input.
'''
class TemplateModelFactory(object):
    def __init__(self, name, description, energies, names_of_parameters,
                 interpolation_degree=1, spline_smoothing_factor=0):
'''
tmf = TemplateModelFactory(mod, mod, ene*uu.keV, ['time'])
tmf.define_parameter_grid('time', time_grid)

# Finally, we loop over our grid and set the interpolation data to the
# template model factory. The units of the fluxes must be a
# differential photon flux!
for ii, tt in enumerate(time_grid):
    tmf.add_interpolation_data(dat[:-1,ii+1], time=tt)

# We can now save our model to disk. The formal is an HDF5 file which
# is saved to the astromodel data directory (~/.astromodels/data). The
# HDF5 file can easily be passed around to other users as all
# information defining the model is stored in the file. The other user
# would place the file in thier astromodels data directory.
tmf.save_data(overwrite=True)
reloaded_table_model = TemplateModel(mod)
reloaded_table_model(ene)
db()
