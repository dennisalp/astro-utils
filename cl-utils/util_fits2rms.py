# Dennis Alp 2017-01-17
# Make an rms map of many .fits images
# On my setup I need to work in the iraf27 environment: source activate iraf27
# Then astropy becomes accesible from both python and python2.7 (both are 2.7.12 anyway)
# time python /Users/silver/Dropbox/bin/util_fits2rms.py

import numpy as np
import os
import pdb

from glob import glob
from pyraf import iraf
from astropy.io import fits
from datetime import date
from scipy.interpolate import griddata
from scipy.ndimage.interpolation import zoom

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

#########
# Parameters
# Logistics
WRK_DIR = "/Users/silver/Dropbox/phd/projects/87a/uplim/rms"
os.chdir(WRK_DIR) #Move to designated directory
files = glob('*.fits') #Find all files.

# Parameters
SNDATE = date(1987, 2, 23)
BOX = [573,599,580,606]
NX = 1150
NY = 1115

#########
# Help functions
# Days for scaling
def get_days(fname): # 1994-09-24_drz_al2.fits
    yr = int(fname[-23:-19])
    month = int(fname[-18:-16].lstrip("0"))
    day = int(fname[-15:-13].lstrip("0"))
    print "Date from filename:", yr, month, day
    return (date(yr, month, day)-SNDATE).days*24*60*60

rms_map = np.zeros((BOX[1]-BOX[0], BOX[3]-BOX[2], len(files)))
counter = 0
dates = np.zeros(len(files))
for img in files[0:]:
    print img[-28:]
    unnormed = fits.open(img)[0].data[BOX[0]:BOX[1],BOX[2]:BOX[3]]
    norm = np.sum(unnormed)
    rms_map[:,:,counter] = unnormed/norm
    
    dates[counter] = get_days(img)
    counter += 1

# Plot detections
rms_map = rms_map.reshape((BOX[1]-BOX[0])*(BOX[3]-BOX[2]),len(files)).T
coefs = np.polyfit(dates,rms_map,1)

for i in range(0,coefs.shape[1]):
    rms_map[:,i] = rms_map[:,i] - np.polyval(coefs[:,i],dates)

rms_map = rms_map.T.reshape(BOX[1]-BOX[0],BOX[3]-BOX[2],len(dates))
rms_map = np.std(rms_map, axis=2)
print np.mean(rms_map)
plt.imshow(rms_map,interpolation='nearest', cmap='afmhot', origin='lower')
plt.show()

rms_map = zoom(rms_map,0.5)
plt.imshow(rms_map,interpolation='nearest', cmap='afmhot', origin='lower')
plt.show()
