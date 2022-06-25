# Dennis Alp 2019-03-15
# Stack all SN 2014J NuSTAR images
# Register images/astrometry
# Print source region files
# Create background light curves
# time python stack.py

import numpy as np
import sys
import os
from pdb import set_trace as db
import time
from glob import glob
from datetime import date
from tqdm import tqdm

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata
from scipy.ndimage.interpolation import rotate
from scipy.ndimage.interpolation import zoom
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import wcs
from astropy import units as au

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)



################################################################
# Help functions
def kev2pi(kev):
    return (kev-1.6)/0.04
def pi2kev(pi):
    return 0.04*pi+1.6

# Just reasonable plot of sky images, good for debugging
def sky_plt(image):
    plt.figure()
    plt.imshow(image,interpolation='nearest', cmap='afmhot', origin='lower')
    plt.show()

# Utility for making fits image out of a numpy array
def mk_fits(image, output, coords):
    hdu = fits.PrimaryHDU(image, header=coords.to_header())
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(output, clobber=True)
    hdulist.close()

def coords2pix(ff, ra, de):
    raref = fits.getval(ff,'TCRVL13', 1)
    radel = fits.getval(ff,'TCDLT13', 1)
    rapix = fits.getval(ff,'TCRPX13', 1)
    deref = fits.getval(ff,'TCRVL14', 1)
    dedel = fits.getval(ff,'TCDLT14', 1)
    depix = fits.getval(ff,'TCRPX14', 1)
    x0 = (ra-raref)*np.cos(np.deg2rad(de))/radel+rapix-1
    y0 = (de-deref)/dedel+depix-1
    return x0, y0

def get_wcs(ff):
    raref = fits.getval(ff,'TCRVL13', 1)
    radel = fits.getval(ff,'TCDLT13', 1)
    rapix = fits.getval(ff,'TCRPX13', 1)
    deref = fits.getval(ff,'TCRVL14', 1)
    dedel = fits.getval(ff,'TCDLT14', 1)
    depix = fits.getval(ff,'TCRPX14', 1)

    coords = wcs.WCS(fits.open(ff)[0].header)
    coords.wcs.crpix = [fov//2, fov//2]
    coords.wcs.cdelt = np.array([scale*radel, scale*dedel])
    coords.wcs.crval = [cc.ra.value, cc.dec.value]
    coords.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    coords.wcs.set_pv([(2, 1, 45.0)])
    
    pixcrd = np.array([[0, 0], [24, 38], [45, 98]], dtype=np.float64)
    world = coords.wcs_pix2world(pixcrd, 0)
    pixcrd2 = coords.wcs_world2pix(world, 0)
    assert np.max(np.abs(pixcrd - pixcrd2)) < 1e-6
    return coords

# This is the combination of the two ULXs in M82, see bachetti14
# "We assume that the persistent emission is composed equally of
# emission from M82 X-1 and X-2 since the contemporaneous Chandra
# imaging shows that these dominate the image and have comparable
# flux."
def get_ulx():
    pulsed = SkyCoord('09h55m51.05s', '69d40m47.90', frame='fk5')
    return SkyCoord(pulsed.ra-4*0.492/3600*au.deg/np.cos(np.deg2rad(pulsed.dec.value)), pulsed.dec+3*0.492/3600*au.deg, frame='fk5')
# We adopt this ('09h55m50.6722s +69d40m49.376s') as the position of the source in NuSTAR

def print_reg(ff, x0, y0, cc):
    out = ff.split('/')[-1]
    out = out.replace('cl.evt', 'src.reg')
    xj = x0-(snj.ra-cc.ra).arcsec/2.46*np.cos(np.deg2rad(snj.dec.value))
    yj = y0+(snj.dec-cc.dec).arcsec/2.46

    regf = open('reg/' + out, 'w') 
    regf.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nphysical\n')
    regf.write('circle(' + str(xj) + ',' + str(yj) + ',' + str(52/2.46) + ')')
    regf.close()



#########
# Parameters
cc = 2.99792458e10 # cm s-1
hh = 6.6260755e-27
DD = 51.2 # kpc
pc = 3.086e18 # cm
kpc = 3.086e21 # cm
mpc = 3.086e24 # cm
kev2erg = 1.60218e-9 # erg keV-1
Msun = 1.989e33 # g
Lsun = 3.828e33 # erg s-1
uu = 1.660539040e-24 # g
only67 = 'only67_' if False else '' # Switch True/False, this helps some string parsing
aligned = '_aligned' if True else '' # Switch True/False, this helps some string parsing

# Logistics
path = '/Users/silver/dat/nus/14j/' + sys.argv[3] + '/'
out_dir = '/Users/silver/box/phd/pro/14j/nus/fig/'
files = sorted(glob(path + '*/nu8???????????01_cl.evt'))
scale = 1
ene_min = kev2pi(float(sys.argv[1])) # keV
ene_max = kev2pi(float(sys.argv[2])) # keV
ns = 1000
snj = SkyCoord('09h55m42.137s', '69d40m25.40s', frame='fk5') # kelly14
# TNS        09:55:42.139 +69:40:26.00
#OSN catalog 09:55:42.12  +69:40:25.9
#OSN catalog 09:55:42.139 +69:40:26.00
#OSN catalog 09:55:42.14  +69:40:26.0

cc = get_ulx()
difference = []

# Allocations
fov = 600
bin_cou = [fov*scale, fov*scale]
bin_lim = [[-fov//2*scale, fov//2*scale], [-fov//2*scale, fov//2*scale]]
img = np.zeros(bin_cou)
bins = np.linspace(0, 1000, ns+1)
spec = np.zeros(ns)
exp_tim = []

imexam = {
    '/Users/silver/dat/nus/14j/' + sys.argv[3] + '/80002092002/nu80002092002A01_cl.evt': (449.69, 474.35),
    '/Users/silver/dat/nus/14j/' + sys.argv[3] + '/80002092002/nu80002092002B01_cl.evt': (448.56, 475.70),
    '/Users/silver/dat/nus/14j/' + sys.argv[3] + '/80002092004/nu80002092004A01_cl.evt': (452.88, 471.07),
    '/Users/silver/dat/nus/14j/' + sys.argv[3] + '/80002092004/nu80002092004B01_cl.evt': (452.14, 472.27),
    '/Users/silver/dat/nus/14j/' + sys.argv[3] + '/80002092006/nu80002092006A01_cl.evt': (446.45, 462.58),
    '/Users/silver/dat/nus/14j/' + sys.argv[3] + '/80002092006/nu80002092006B01_cl.evt': (445.95, 463.88),
    '/Users/silver/dat/nus/14j/' + sys.argv[3] + '/80002092007/nu80002092007A01_cl.evt': (461.68, 489.87),
    '/Users/silver/dat/nus/14j/' + sys.argv[3] + '/80002092007/nu80002092007B01_cl.evt': (461.20, 491.64),
    '/Users/silver/dat/nus/14j/' + sys.argv[3] + '/80002092008/nu80002092008A01_cl.evt': (459.36, 495.89),
    '/Users/silver/dat/nus/14j/' + sys.argv[3] + '/80002092008/nu80002092008B01_cl.evt': (459.27, 497.55),
    '/Users/silver/dat/nus/14j/' + sys.argv[3] + '/80002092009/nu80002092009A01_cl.evt': (447.96, 488.33),
    '/Users/silver/dat/nus/14j/' + sys.argv[3] + '/80002092009/nu80002092009B01_cl.evt': (447.61, 489.85),
    '/Users/silver/dat/nus/14j/' + sys.argv[3] + '/80002092011/nu80002092011A01_cl.evt': (425.22, 500.43),
    '/Users/silver/dat/nus/14j/' + sys.argv[3] + '/80002092011/nu80002092011B01_cl.evt': (425.71, 501.96)}

################################################################
# Main
for ii, ff in enumerate(tqdm(files)):
    if only67 == 'only67_' and not '80002092006' in ff and not '80002092007' in ff:
        continue
    
    x1, y1 = coords2pix(ff, cc.ra.degree, cc.dec.degree)
    x0, y0 = imexam[ff]
    print_reg(ff, x0, y0, cc)
    difference.append([x0-x1, y0-y1])
    if aligned == '':
        x0, y0 = coords2pix(ff, cc.ra.degree, cc.dec.degree)
        
    dat = fits.open(ff)[1].data
    tt = dat['TIME']
    ene = dat['PI']
    xx = dat['X']
    yy = dat['Y']

    xx = xx-x0
    yy = yy-y0
    exp_tim.append(fits.getval(ff, 'EXPOSURE', 1))
    
    # Make Image
    hig_ene = (ene > ene_min) & (ene < ene_max)
    tmp, xedges, yedges = np.histogram2d(scale*xx[hig_ene], scale*yy[hig_ene], bins=bin_cou, range=bin_lim)
    img = img + tmp

    # Make spectrum
    idx = (np.abs(xx) < fov//2) & (np.abs(yy) < fov//2) & hig_ene
    tmp, _ = np.histogram(pi2kev(ene[idx]), bins)
    spec = spec + tmp

    # Make background light curves
    fig = plt.figure(figsize=(15, 3.75))
    bkg = (ene > ene_min) & (ene < ene_max) & (xx**2+yy**2 > 30**2)
    tbin = int(np.ceil((tt[bkg].max()-tt[bkg].min())/256.))
    plt.hist(tt[bkg], tbin, histtype='step')
    fig.savefig('fig/' + ff.split('/')[-1].replace('cl.evt', 'bkg_lc_') + sys.argv[1] + '-' + sys.argv[2] + '_' + path.split('/')[-2] + '.pdf',bbox_inches='tight', pad_inches=0.01)
#    plt.show()
    plt.close()



################################################################
# Spectrum
tmp = np.zeros(2*ns+2)
tmp[::2] = bins
tmp[1::2] = bins
bins = tmp.copy()
tmp = np.zeros(2*ns+2)
tmp[1:-2:2] = spec
tmp[2:-1:2] = spec
tmp[0] = tmp[-1] = 0
spec = tmp.copy()
plt.loglog(bins, spec, 'k')

# Image
my_wcs = get_wcs(files[0])
mk_fits(img.T, out_dir + 'stack_' + only67 + sys.argv[1] + '-' + sys.argv[2] + aligned + '.fits', my_wcs)
print('Exposure time:', np.sum(np.array(exp_tim)), 'Individual:', np.array(exp_tim))
print('Counts:', np.sum(img))
#sky_plt(img.T)

difference = np.array(difference)
db()

#  3-79 keV PEAK  FWHM
#default 2515.26  8.57
#aligned 2652.34  8.19
# 20-79 keV PEAK MFWHM
#default   44.88 11.18
#aligned   47.20 10.12
