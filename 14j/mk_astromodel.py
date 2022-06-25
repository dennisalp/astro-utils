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
# First, we grab a function and make an energy grid.
model = Band()
# we won't need to modify the normalization
model.K = 1.

# if no units are provided for the energy grid, keV will be assumed!
energies = np.logspace(1, 3, 50)

# Now we define a template model factory. This takes a name, a
# description, the energy grid and an array of parameter names as
# input.

tmf = TemplateModelFactory('my_template', 'A test template', energies, ['alpha', 'xp', 'beta'])
# Now, define our grid in parameter space. While we are using a
# function here, this grid could be from a text file, a database of
# simulations, etc. We then assign these grid points to the template
# model factory.

alpha_grid = np.linspace(-1.5, 1, 15)
beta_grid = np.linspace(-3.5, -1.6, 15)
xp_grid = np.logspace(1, 3, 20)

tmf.define_parameter_grid('alpha', alpha_grid)
tmf.define_parameter_grid('beta', beta_grid)
tmf.define_parameter_grid('xp', xp_grid)

# Finally, we loop over our grid and set the interpolation data to the
# template model factory. The units of the fluxes must be a
# differential photon flux!
for a in alpha_grid:
    for b in beta_grid:
        for xp in xp_grid:
            # change our model parameters
            model.alpha = a
            model.beta = b
            model.xp = xp
            tmf.add_interpolation_data(model(energies), alpha=a, xp=xp, beta=b)



################################################################
# We can now save our model to disk. The formal is an HDF5 file which
# is saved to the astromodel data directory (~/.astromodels/data). The
# HDF5 file can easily be passed around to other users as all
# information defining the model is stored in the file. The other user
# would place the file in thier astromodels data directory.
tmf.save_data(overwrite=True)
reloaded_table_model = TemplateModel('my_template')
reloaded_table_model(energies)
