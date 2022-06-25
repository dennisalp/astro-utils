# Dennis Alp 2016-12-11
# This is supposed to merge colors of many observations
# time python /Users/$USER/Dropbox/bin/util_fits2col.py

import numpy as np
import os
import pdb

from glob import glob
from astropy.io import fits

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

#########
# Parameters
WRK_DIR = "/Users/silver/Dropbox/phd/projects/87a/fold/"
BOX=[940,1620,2988,3668]
COLGRAD = 0.8 #1 white #1.25 for red
ARCGRAD = 300 #1e3 white # 4000 for red
REDLIN = 70.
REDPOW = 0.8
CAP = 0.06
DPI = 870
OUTPUT = "white.png"
INTERPOLATION = 'spline36'
INTERPOLATION = 'nearest'
methods = [None, 'none', 'nearest', 'bilinear', 'bicubic', 'spline16',
           'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric',
           'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']

os.chdir(WRK_DIR) #Move to designated directory
files = glob('*.fits') #Find all files.
NN = len(files)
output = np.zeros((1115,1150,3))
#output = np.zeros((550,550,3))

colors = {
"alm212014-09-02_sum_al2.fits": [255,255,255],
"alm542014-09-02_sum_al2.fits": [255,255,255],
"alm652014-09-02_sum_al2.fits": [255,255,255],
"w32252009-12-13_drz_al2.fits": [148,  0,212],
"w33362009-12-13_drz_al2.fits": [148,  0,255],
"w35022016-06-08_drz_al2.fits": [  0,177,118],
"w35552009-12-13_drz_al2.fits": [  0,255,  0],
"w36452016-06-08_drz_al2.fits": [255,  0,  0],
"w38142009-12-13_drz_al2.fits": [127,  0,  0],
"w3_b_2015-05-24_drz_al2.fits": [100,  0,255],
"w3_b_2016-06-08_drz_al2.fits": [100,  0,255],
"w3_r_2015-05-24_drz_al2.fits": [255,103,  0],
"w3_r_2016-06-08_drz_al2.fits": [255,103,  0],
}

weights = {
"alm212014-09-02_sum_al2.fits": 0,
"alm542014-09-02_sum_al2.fits": 0,
"alm652014-09-02_sum_al2.fits": 0,
"w32252009-12-13_drz_al2.fits": 1, # violet
"w33362009-12-13_drz_al2.fits": 1, # violet
"w35022016-06-08_drz_al2.fits": 0, # null
"w35552009-12-13_drz_al2.fits": 4, # green
"w36452016-06-08_drz_al2.fits": 0, # null
"w38142009-12-13_drz_al2.fits": 1, # dark red
"w3_b_2015-05-24_drz_al2.fits": 3, # blue
"w3_b_2016-06-08_drz_al2.fits": 3, # blue
"w3_r_2015-05-24_drz_al2.fits": 1.75, # red
"w3_r_2016-06-08_drz_al2.fits": 1.75, # red
}

#white
colors = {
"alm212014-09-02_sum_al2.fits": [255,255,255],
"alm542014-09-02_sum_al2.fits": [255,255,255],
"alm652014-09-02_sum_al2.fits": [255,255,255],
"w32252009-12-13_drz_al2.fits": [  0,  0,255],
"w33362009-12-13_drz_al2.fits": [  0,  0,255],
"w35022016-06-08_drz_al2.fits": [  0,177,118],
"w35552009-12-13_drz_al2.fits": [  0,255,  0],
"w36452016-06-08_drz_al2.fits": [255,  0,  0],
"w38142009-12-13_drz_al2.fits": [255,  0,  0],
"w3_b_2015-05-24_drz_al2.fits": [  0, 64,255],
"w3_b_2016-06-08_drz_al2.fits": [  0, 64,255],
"w3_r_2015-05-24_drz_al2.fits": [255,  0,  0],
"w3_r_2016-06-08_drz_al2.fits": [255,  0,  0],
}

weights = {
"alm212014-09-02_sum_al2.fits": 0,
"alm542014-09-02_sum_al2.fits": 0,
"alm652014-09-02_sum_al2.fits": 0,
"w32252009-12-13_drz_al2.fits": 2.3, # violet
"w33362009-12-13_drz_al2.fits": 2.3, # violet
"w35022016-06-08_drz_al2.fits": 0, # null
"w35552009-12-13_drz_al2.fits": 2.4, # green 4
"w36452016-06-08_drz_al2.fits": 0, # null
"w38142009-12-13_drz_al2.fits": 1, # dark red
"w3_b_2015-05-24_drz_al2.fits": 2, # blue 3
"w3_b_2016-06-08_drz_al2.fits": 2, # blue
"w3_r_2015-05-24_drz_al2.fits": 1, # red 1.75
"w3_r_2016-06-08_drz_al2.fits": 1, # red
}

    
#red
#colors = {
#"alm212014-09-02_sum_al2.fits": [255,255,255],
#"alm542014-09-02_sum_al2.fits": [255,255,255],
#"alm652014-09-02_sum_al2.fits": [255,255,255],
#"w32252009-12-13_drz_al2.fits": [  0,  0,212],
#"w33362009-12-13_drz_al2.fits": [  0,  0,255],
#"w35022016-06-08_drz_al2.fits": [  0,  0,  0],
#"w35552009-12-13_drz_al2.fits": [  0,255,  0],
#"w36452016-06-08_drz_al2.fits": [  0,  0,  0],
#"w38142009-12-13_drz_al2.fits": [127,0,  0],
#"w3_b_2015-05-24_drz_al2.fits": [ 0, 0,255],
#"w3_b_2016-06-08_drz_al2.fits": [ 0, 0,255],
#"w3_r_2015-05-24_drz_al2.fits": [255,0,  0],
#"w3_r_2016-06-08_drz_al2.fits": [255,0,  0],
#}
#    
#weights = {
#"alm212014-09-02_sum_al2.fits": 0,
#"alm542014-09-02_sum_al2.fits": 0,
#"alm652014-09-02_sum_al2.fits": 0,
#"w32252009-12-13_drz_al2.fits": 1, # violet
#"w33362009-12-13_drz_al2.fits": 1, # violet
#"w35022016-06-08_drz_al2.fits": 0, # null
#"w35552009-12-13_drz_al2.fits": 6, # green
#"w36452016-06-08_drz_al2.fits": 0, # null
#"w38142009-12-13_drz_al2.fits": 1, # dark red
#"w3_b_2015-05-24_drz_al2.fits": 1, # blue
#"w3_b_2016-06-08_drz_al2.fits": 1, # blue
#"w3_r_2015-05-24_drz_al2.fits": 7, # red
#"w3_r_2016-06-08_drz_al2.fits": 7, # red
#}

#blue
#colors = {
#"alm212014-09-02_sum_al2.fits": [255,255,255],
#"alm542014-09-02_sum_al2.fits": [255,255,255],
#"alm652014-09-02_sum_al2.fits": [255,255,255],
#"w32252009-12-13_drz_al2.fits": [  14,  70,255],
#"w33362009-12-13_drz_al2.fits": [  14,  70,255],
#"w35022016-06-08_drz_al2.fits": [  14,  70,255],
#"w35552009-12-13_drz_al2.fits": [  14,  70,255],
#"w36452016-06-08_drz_al2.fits": [  14,  70,255],
#"w38142009-12-13_drz_al2.fits": [  14,  70,255],
#"w3_b_2015-05-24_drz_al2.fits": [  14,  70,255],
#"w3_b_2016-06-08_drz_al2.fits": [  14,  70,255],
#"w3_r_2015-05-24_drz_al2.fits": [  14,  70,255],
#"w3_r_2016-06-08_drz_al2.fits": [  14,  70,255],
#}
#    
#weights = {
#"alm212014-09-02_sum_al2.fits": 0,
#"alm542014-09-02_sum_al2.fits": 0,
#"alm652014-09-02_sum_al2.fits": 0,
#"w32252009-12-13_drz_al2.fits": 1, # violet
#"w33362009-12-13_drz_al2.fits": 1, # violet
#"w35022016-06-08_drz_al2.fits": 0, # null
#"w35552009-12-13_drz_al2.fits": 1, # green
#"w36452016-06-08_drz_al2.fits": 0, # null
#"w38142009-12-13_drz_al2.fits": 1, # dark red
#"w3_b_2015-05-24_drz_al2.fits": 1, # blue
#"w3_b_2016-06-08_drz_al2.fits": 1, # blue
#"w3_r_2015-05-24_drz_al2.fits": 1, # red
#"w3_r_2016-06-08_drz_al2.fits": 1, # red
#}

#########
# Main fold
for ff in files[0:]:
    key = ff[-28:]

    # Do main part
    ff = fits.open(ff)[0].data
#    ff = fits.open(ff)[0].data[BOX[0]:BOX[2],BOX[1]:BOX[3]]
        
    goodmin = np.min(ff[450:650,450:650]) # Find a good minimum value
    ff = np.where(ff < goodmin, goodmin, ff) # Set lower bound
    ff = ff - np.min(ff) # Move lower bound to zero
    ff = np.where(ff > 120., 120., ff) # Set lower bound to WFC3 saturation
    norm = np.sum(ff[350:600,700:850]) # Find normalization
    ff = ff/norm
    print(key,norm)
    
    # Add the layers
    for col in range(0,3):
        output[:,:,col] = output[:,:,col] + weights[key]*ff*colors[key][col]

    # Save some for reddenning
    if key == "w3_b_2015-05-24_drz_al2.fits":
        b15 = ff
    if key == "w3_b_2016-06-08_drz_al2.fits":
        b16 = ff
    elif key == "w3_r_2015-05-24_drz_al2.fits":
        r15 = ff
    elif key == "w3_r_2016-06-08_drz_al2.fits":
        r16 = ff

#########
# Reddenning 
#redness = (output[:,:,0]/output[:,:,2])
#output[:,:,0] = np.where(redness > 1., redness**REDDEN*output[:,:,0], output[:,:,0])
delta = r15+r16-b15-b16
#norm = np.sum(delta[350:600,700:850]) # Find normalization
#delta = delta/norm
delta = np.where(delta > 0, delta, 0.)
output[:,:,0] = output[:,:,0] + REDLIN*delta**REDPOW

#########
# Post production
output = output/np.max(output)

# Reduce background
#output = np.where(output < 1e-4, 12, output)

# Colorscale
output = np.arctan(output**COLGRAD*ARCGRAD)/np.pi*2

# Plotting
plt.imshow(output, interpolation=INTERPOLATION, origin='lower')
#plt.savefig(OUTPUT, dpi=DPI, bbox_inches='tight', pad_inches=0)
plt.show()

# Crop
#crop = plt.imread(OUTPUT)
#plt.imsave("crop_" + OUTPUT,crop[BOX[0]:BOX[2],BOX[1]:BOX[3],:])
#crop[BOX[0]:BOX[2],BOX[1]:BOX[3],:] = np.zeros((2048,2048,4))
#plt.imsave("check_" + OUTPUT,crop)
