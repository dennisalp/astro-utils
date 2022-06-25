'''
2017-07-01, Dennis Alp, dalp@kth.se

Compute the area of the source and background spectra extraction regions that maximize the S/N. This requires the PSF/encircled energy and is a function of the source and background areas. See ../art/main.pdf. This script also extract the effective area from an .arf (scaled to no PSF correction).
'''

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
from scipy.integrate import cumtrapz
from scipy.ndimage.interpolation import rotate
from scipy.ndimage.interpolation import zoom
from scipy.optimize import curve_fit
from astropy.io import fits
from astropy.coordinates import SkyCoord

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
def mk_fits(image, output):
    hdu = fits.PrimaryHDU(image)
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

# Logistics
path = '/Users/silver/dat/nus/14j/v01/'
out_dir = '/Users/silver/box/phd/pro/14j/nus/fig/'
psf1 = np.loadtxt('/Users/silver/box/sci/lib/m/madsen15b/fig15.txt')
psf2 = np.loadtxt('/Users/silver/box/sci/lib/m/madsen15b/nustar_observatory_guide_1p0.txt')
files = sorted(glob(path + '*/nu????????????01_cl.evt'))
scale = 1
ene_min = kev2pi(10.) # keV
ene_max = kev2pi(78.) # keV
ns = 1000
cc = SkyCoord('09h55m42.217s', '69d40m26.56s', frame='fk5')

# Allocations
fov = 600
bin_cou = [fov*scale, fov*scale]
bin_lim = [[-fov//2*scale, fov//2*scale], [-fov//2*scale, fov//2*scale]]
img = np.zeros(bin_cou)
exp_tim = 0



################################################################
# Main
for ii, ff in enumerate(tqdm(files)):
    x0, y0 = coords2pix(ff, cc.ra.degree, cc.dec.degree)
    dat = fits.open(ff)[1].data
    ene = dat['PI']
    xx = dat['X']
    yy = dat['Y']

    xx = xx-x0
    yy = yy-y0
    exp_tim += fits.getval(ff, 'EXPOSURE', 1)
    
    # Make Image
    hig_ene = (ene > ene_min) & (ene < ene_max)
    tmp, xedges, yedges = np.histogram2d(scale*xx[hig_ene], scale*yy[hig_ene], bins=bin_cou, range=bin_lim)
    img = img + tmp


    
################################################################
# Fit PSF
print('Exposure time:', exp_tim)
print('Counts:', np.sum(img))
#sky_plt(img.T)
img = np.where(img == 0, 6, img)

def gauss(dummy, x0, y0, a1, s1, a2, s2, a3, s3, bb):
    print('{0:8.2f}{1:8.2f}{2:8.2f}{3:8.2f}{4:8.2f}{5:8.2f}{6:8.2f}{7:8.2f}{8:8.2f}'.format(x0, y0, a1, s1, a2, s2, a3, s3, bb))
    xx = np.arange(0,fov)
    yy = np.arange(0,fov)[:,np.newaxis]
    xx = xx-x0
    yy = yy-y0
    rr = xx**2+yy**2
    global fit_gau
    fit_gau = a1*np.exp(-0.5*rr/s1**2)+a2*np.exp(-0.5*rr/s2**2)+a3*np.exp(-0.5*rr/s3**2)+bb
    return np.reshape(fit_gau, fov**2)

guess = np.array([fov//2,    fov//2,   200,   4,   20,   8,   300,   2,   6])
bl = np.array([fov//2-20, fov//2-20, 200/3, 4/2, 20/3, 8/3, 300/3, 2/3, 6/3])
bu = np.array([fov//2+20, fov//2+20, 200*3, 4*2, 20*3, 8*3, 300*3, 2*3, 6*3])
pp, covar = curve_fit(gauss, np.arange(0, fov**2), np.reshape(img, fov**2), guess, bounds=(bl, bu))

res = img-fit_gau
rr = np.arange(1,500,0.1)
radial_profile = pp[2]*np.exp(-rr**2/(2*(pp[3]*2.46)**2))+pp[4]*np.exp(-rr**2/(2*(pp[5]*2.46)**2))+pp[6]*np.exp(-rr**2/(2*(pp[7]*2.46)**2))
psf3 = cumtrapz(radial_profile*rr, rr, initial=0.)
psf3 = psf3/psf3[-1]

fig = plt.figure(figsize=(5, 3.75))
plt.semilogx(psf1[:,0], psf1[:,1])
plt.semilogx(psf2[:,0], psf2[:,1])
plt.semilogx(rr, psf3)
plt.xlabel('Radius (arcsec)')
plt.ylabel('Encircled energy')
plt.legend(['madsen15b','NuSTAR Observatory Guide 1.0', 'My fit to M82 (not even a single point source), see bachetti14'])
fig.savefig('fig/psf_frac.pdf',bbox_inches='tight', pad_inches=0.01)

header = 'Encircled energy at roughly 10 keV (almost constant at higher energies) for NuSTAR on-axis. Obtained by fits to the TWO (almost co-spatial) sources in M82. The first column is the radius in arcsec and the second column is the encircled energy.'
np.savetxt('nus/my_psf_frac.txt', np.c_[rr,psf3], header=header)
header = 'Encircled energy at 10 keV (almost constant at higher energies) for NuSTAR on-axis. see Fig 15 of madsen15b. The first column is the radius in arcsec and the second column is the encircled energy.'
np.savetxt('nus/madsen15b.txt', psf1, header=header)
header = 'Encircled energy at 10 keV (almost constant at higher energies) for NuSTAR on-axis. see NuSTAR Observatory Guide 1.0. The first column is the radius in arcsec and the second column is the encircled energy.'
np.savetxt('nus/nustar_observatory_guide_1p0.txt', psf2, header=header)
#plt.show()

# Save files and get the arf
fig = plt.figure(figsize=(5, 3.75))
arf = fits.open('/Users/silver/dat/nus/14j/v03/80002092007/nu80002092007A01_sr.arf')[1].data
arf = np.c_[arf.field(0)[:,np.newaxis], arf.field(2)[:,np.newaxis]/griddata(psf1[:,0], psf1[:,1], 49.2)]
plt.loglog(arf[:,0], arf[:,1])
fig.savefig('fig/effective_area.pdf',bbox_inches='tight', pad_inches=0.01)
header = 'Effective area (cm^2, second column) for one of the telescopes of NuSTAR. Not corrected for contained PSF fraction/aperture correction. The first column is the energy at the lower end of the bin in keV.'
np.savetxt('nus/effective_area.txt', arf, header=header)
#plt.show()



################################################################
# Compute S/N
def get_psf(rr):
    return griddata(psf1[:,0], psf1[:,1], rr, method='cubic', fill_value=np.nan)

# work in radius instead
# v01, v02, and v03 are using src rad 50" and bkg rad 84"
src_rad = np.arange( 2, 100, 0.5)[:,np.newaxis]
bkg_rad = np.arange(50, 200, 1)[np.newaxis,:]
src_rad, bkg_rad = np.meshgrid(src_rad, bkg_rad, indexing='ij')
psf_frac = get_psf(src_rad)
surf = psf_frac/np.sqrt(src_rad+src_rad**2/bkg_rad)
sri, bki = divmod(surf.argmax(), surf.shape[1])

fig = plt.figure(figsize=(5, 3.75))
plt.pcolor(src_rad, bkg_rad, surf)
plt.plot(src_rad[sri,bki], bkg_rad[sri,bki], 'o', color='gray')
print('(arcsec)')
print('For arbitrary total area:{0:4.0f} src rad, {1:4.0f} bkg rad.'.format(src_rad[sri,bki], bkg_rad[sri,bki]))

plt.plot(src_rad[:,0], np.sqrt(330**2/np.pi-src_rad[:,0]**2), 'w') # Roughly the available area
line = psf_frac[:,0]/np.sqrt(src_rad[:,0]+src_rad[:,0]**2/np.sqrt(330**2/np.pi-src_rad[:,0]**2))
ii = line.argmax()
plt.plot(src_rad[ii,0], np.sqrt(330**2/np.pi-src_rad[ii,0]**2), 'wo')
print('Constrained to detector area:{0:4.0f} src rad, {1:4.0f} bkg rad-equivalent.'.format(src_rad[ii,0], np.sqrt(330**2/np.pi-src_rad[ii,0]**2)))

plt.xlabel('Source radius (arcsec)')
plt.ylabel('Background radius-equivalent (arcsec)')
plt.colorbar()
fig.savefig('fig/src_bkg_reg.pdf',bbox_inches='tight', pad_inches=0.01)
plt.show()
db()
