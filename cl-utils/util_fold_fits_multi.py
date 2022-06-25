# Dennis Alp 2016-03-01
# This is supposed to fold many fits files into an average and takes care of some headers
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
WRK_DIR = "/Users/silver/Dropbox/phd/data/alm/87a/all"
os.chdir(WRK_DIR) #Move to designated directory
files = glob('*nat.image*') #Find all files.
files.append('c23.icrs.cont213.b2.multi7.fits')
OUTPUT = "alma_fold.fits"
NN = len(files)

# Allocate
folded = np.zeros(fits.open(files[0])[0].data[0].shape)
selected = [[0,1,59,60,61,62,63],[0,1,2,3,4,5,6,7,8,56,57,58,59,60,61,62,63]]
for img_path in files:
    dat = fits.open(img_path)[0].data
    if 'co21.c23' in img_path:
        folded += np.mean(dat[selected[0],:,:], axis=0)
        for i in selected[0]:
            plt.imshow(dat[i,270:330,300:360], interpolation='nearest', cmap='afmhot', origin='lower')
#            plt.imshow(dat[i], interpolation='nearest', cmap='afmhot', origin='lower')
            plt.plot(32,21,'bo')
            plt.show()
    if 'sio54.c23' in img_path:
        folded += np.mean(dat[selected[1],:,:], axis=0)
        for i in selected[1]:
            plt.imshow(dat[i,270:330,300:360], interpolation='nearest', cmap='afmhot', origin='lower')
#            plt.imshow(dat[i], interpolation='nearest', cmap='afmhot', origin='lower')
            plt.plot(32,21,'bo')
            plt.show()
    if 'c23.icrs' in img_path:
        folded += dat[0,0]

hdu = fits.PrimaryHDU(folded)
hdulist = fits.HDUList([hdu])
hdulist.writeto(OUTPUT, clobber=True)
