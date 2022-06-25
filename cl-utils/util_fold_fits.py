# Dennis Alp 2016-12-11
# This is supposed to fold many fits files into an average and takes care of some headers
# It assumes unit to be electrons (counts) per second
# time python /Users/$USER/Dropbox/bin/util_fold_fits.py

import numpy as np
import os
import pdb

from glob import glob
from astropy.io import fits

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

#########
# Parameters
WRK_DIR = "/Users/silver/Dropbox/phd/projects/87a/uplim/data/fold/"
os.chdir(WRK_DIR) #Move to designated directory
files = glob('*.fits') #Find all files.
NN = len(files)
INTERPOLATION = 'nearest'
DATE_OBS = "0306-00-00"
OUTPUT = "../ac_r_" + DATE_OBS + "_drz_al2.fits"

# Allocate
folded = np.zeros(fits.open(files[0])[0].shape)
exposure = 0

for img_path in files[0:]:
    dat = fits.open(img_path)[0].data
    duration = fits.getval(img_path,'EXPTIME')
    folded = folded + dat*duration
    exposure += duration
    print duration

folded = folded/exposure
hdu = fits.PrimaryHDU(folded)
hdulist = fits.HDUList([hdu])
prihdr = hdulist[0].header
prihdr['EXPTIME'] = exposure
prihdr['DATE-OBS'] = DATE_OBS
hdulist.writeto(OUTPUT, clobber=True)
