'''
2019-06-15, Dennis Alp, dalp@kth.se

Find the position of SN 1987A in HST/FOC images from 1990 and register to Gaia DR2.

Use the IRAF task ccmap to register the images onto the Gaia solution. 
ccmap @all.inp all.dat solutions=@all.sol images=@all.fits results=@all.res update=yes insystem=icrs lngunits=degrees latunits=degrees interactive=no
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LogNorm
import os
import time
from pdb import set_trace as db

from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units
from scipy.optimize import curve_fit
from scipy.interpolate import griddata
from scipy.stats import sem
from glob import glob

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)



################################################################
################################################################



def gen_mask(x0, y0, sig):
    xx = np.arange(0,nx)
    yy = np.arange(0,ny)[:,none]
    xx = xx-x0
    yy = yy-y0
    return np.sqrt(xx**2+yy**2) < 1.85*sig
 
def gauss(dummy, x0, y0, a1, s1, a2, s2, bb):
    xx = np.arange(0,nx)
    yy = np.arange(0,ny)[:,np.newaxis]
    xx = xx-x0
    yy = yy-y0
    rr = xx**2+yy**2
    global fit_gau
    fit_gau = a1*np.exp(-0.5*rr/s1**2)+a2*np.exp(-0.5*rr/s2**2)+bb
    return np.reshape(fit_gau, nx*ny)

def pix2coords(ff, xx, yy):
    rapix = fits.getval(ff,'CRPIX1', 0)
    depix = fits.getval(ff,'CRPIX2', 0)
    xx = xx - rapix
    yy = yy - depix

    cd11 = fits.getval(ff,'CD1_1', 0)
    cd12 = fits.getval(ff,'CD1_2', 0)
    cd21 = fits.getval(ff,'CD2_1', 0)
    cd22 = fits.getval(ff,'CD2_2', 0)
    ra = cd11*xx+cd12*yy
    de = cd21*xx+cd22*yy

    raref = fits.getval(ff,'CRVAL1', 0)
    deref = fits.getval(ff,'CRVAL2', 0)
    de = de+deref
    ra = ra/np.cos(np.deg2rad(de))+raref    
    return ra, de 



################################################################
################################################################



#########
# Parameters
dat_dir = "dat/"
interpolation = 'nearest'
sna = SkyCoord('5h35m27.9875s','-69d16m11.107s', frame='icrs')

files = ['x0c80102t_c1f.fits',
         'x0c80103t_c1f.fits',
         'x0c80104t_c1f.fits',
         'x0c80105t_c1f.fits',
         'x0c80106t_c1f.fits',
         'x0c80107t_c1f.fits']

nn = len(files)
cc = []
pars = []

nx = 31
ny = 31

#########
# Gaussian parameters
x0 = np.array([ 282., 281., 281., 276., 293., 292.])-1 # Subtract 1 to map from image to Python
y0 = np.array([ 251., 251., 251., 256., 251., 251.])-1
a1 = np.array([  40.,  40., 300., 400.,  36.,  36.])
s1 = np.array([   2.,   2.,   2.,   2.,   2.,   2.])
a2 = np.array([  10.,  10.,  80.,  80.,   8.,   8.])
s2 = np.array([   4.,   4.,   4.,   4.,   4.,   4.])
bb = np.array([   3.,   3.,  15.,  20.,   4.,   4.])

#########
# Fit Gaussians for all observations
for ii, img_path in enumerate(files):
# Gauss
    img = fits.open(dat_dir + img_path)[0].data
    img = img[int(y0[ii])-ny//2:int(y0[ii])+ny//2+1, int(x0[ii])-nx//2:int(x0[ii])+nx//2+1]
 
    guess = np.array([nx//2, ny//2, a1[ii], s1[ii], a2[ii], s2[ii], bb[ii]])
    bl = np.array([nx//2-5, ny//2-5, a1[ii]/3, s1[ii]/3, a2[ii]/3, s2[ii]/3, bb[ii]/3])
    bu = np.array([nx//2+5, ny//2+5, a1[ii]*3, s1[ii]*3, a2[ii]*3, s2[ii]*3, bb[ii]*3])
    pp, covar = curve_fit(gauss, np.arange(0, nx*ny), np.reshape(img, nx*ny), guess, bounds=(bl, bu))

    pp[0] = pp[0]+x0[ii]-nx//2+1 # Add 1 to map from Python to image
    pp[1] = pp[1]+y0[ii]-ny//2+1
    print("\n" + img_path + "\nxx = {0:17.13f}\nyy = {1:17.13f}\na1 = {2:17.13f}\ns1 = {3:17.13f}\na2 = {4:17.13f}\ns2 = {5:17.13f}\nbb = {6:17.13f}".format(pp[0], pp[1], pp[2], pp[3], pp[4], pp[5], pp[6]))
    pars.append(pp)
 
# Plots
    fig = plt.figure(figsize=(5, 3.75))
    plt.imshow(img, cmap='afmhot', origin='lower', interpolation=interpolation)
    plt.plot(pp[0]-1-x0[ii]+nx//2, pp[1]-1-y0[ii]+ny//2, 'ko')
    fig.savefig('plt/' + img_path + '_img.pdf',bbox_inches='tight', pad_inches=0.01)
#    plt.show()
    plt.close(fig)
 
    fig = plt.figure(figsize=(5, 3.75))
    plt.imshow(img-fit_gau, cmap='afmhot', origin='lower', interpolation=interpolation)
    fig.savefig('plt/' + img_path + '_res.pdf',bbox_inches='tight', pad_inches=0.01)
#    plt.show()
    plt.close(fig)

    ra, de = pix2coords(dat_dir+img_path, pp[0], pp[1])
    cc.append(np.array([ra, de]))


    
################################################################
################################################################



cc = np.array(cc)
c0 = cc[:,0].mean()
c1 = cc[:,1].mean()
pars = np.array(pars)
std = np.std(cc,axis=0)
ra2mas = np.cos(np.deg2rad(cc[:,1].mean()))*3600e3
std[0] = std[0]*ra2mas
std[1] = std[1]*3600e3
sna = SkyCoord('5h35m27.9875s','-69d16m11.107s', frame='icrs')
sn0, sn1 = sna.ra.value, sna.dec.value
ra = np.average((c0, sn0))
de = np.average((c1, sn1))
new = SkyCoord(ra*units.degree, de*units.degree, frame='icrs')
print('\nThe position of SN 1987A is ' + new.to_string('hmsdms') + '.')
tmp = np.round(20/2**0.5/np.cos(np.deg2rad(de)), 2)
print('We assume a 1-sigma uncertainty of 14 mas in both RA and dec (' + str(tmp) + '\" in RA), resulting in an aboslute uncertainty of 20 mas.\n')

col = ['purple', 'purple', 'b', 'g', 'gold', 'gold']
fmt = ['o', 's', 'o', 'o', 'o', 's']
fig = plt.figure(figsize=(5, 3.75))
plt.imshow(img, cmap='afmhot', origin='lower', interpolation=interpolation)
for ii, pp in enumerate(pars):
    plt.plot(pp[0]-1-x0[ii]+nx//2, pp[1]-1-y0[ii]+ny//2, color=col[ii], marker=fmt[ii])

fig.savefig('plt/map.pdf',bbox_inches='tight', pad_inches=0.01)
plt.close(fig)

fig = plt.figure(figsize=(5, 3.75))
ax = plt.gca()
for ii, c in enumerate(cc):
    plt.plot((c[0]-ra)*ra2mas, (c[1]-de)*3600e3, color=col[ii], marker=fmt[ii])

plt.plot((sn0-ra)*ra2mas, (sn1-de)*3600e3, 'x', color='gray', ms=10, mew=1.5)
plt.plot((c0-ra)*ra2mas, (c1-de)*3600e3, 'x', color='crimson', ms=10, mew=1.5)
plt.plot(0, 0, 'x', color='k', ms=10, mew=3)
ellipse = patches.Ellipse(xy=((c0-ra)*ra2mas, (c1-de)*3600e3), width=2*std[0], height=2*std[1], edgecolor='crimson', fc='None', lw=1.5)
ax.add_artist(ellipse)
ellipse = patches.Ellipse(xy=((sn0-ra)*ra2mas, (sn1-de)*3600e3), width=2*6, height=2*4, edgecolor='gray', fc='None', lw=1.5)
ax.add_artist(ellipse)
ellipse = patches.Ellipse(xy=(0, 0), width=2*20/2**0.5, height=2*20/2**0.5, edgecolor='k', fc='None', lw=3)
ax.add_artist(ellipse)
#ellipse.set_zorder(5)
ax.set_xlim([-16,16])
ax.set_ylim([-16,16])
ax.invert_xaxis()
ax.set_aspect('equal')
plt.yticks([-10, 0, 10])
plt.ylabel('$\Delta\delta$ (mas)')
plt.xlabel('$\Delta\\alpha$ (mas)')
fig.savefig('plt/sky.pdf',bbox_inches='tight', pad_inches=0.01)
plt.show()
db()
plt.close(fig)

#x0c80102t_c1f.fits
#iraf: 281.54, 250.86
#xx = 281.6819540685120
#yy = 250.8372337828310
#a1 =  32.9092908771487
#s1 =   2.4881306301681
#a2 =   8.6519914487003
#s2 =   5.6580090642917
#bb =   4.9297163091995
#
#x0c80103t_c1f.fits
#iraf: 280.76, 251.30
#xx = 280.8713630769697
#yy = 251.3987576039622
#a1 =  30.7527072594549
#s1 =   2.4095941437364
#a2 =   9.7056104905877
#s2 =   5.6158639571567
#bb =   5.4226913269846
#
#x0c80104t_c1f.fits
#iraf: 280.33, 251.50
#xx = 280.5088098702037
#yy = 251.5723671714402
#a1 = 288.0343505977374
#s1 =   2.8015862792391
#a2 =  73.4816141960273
#s2 =   6.0506308839302
#bb =  22.4670836314193
#
#x0c80105t_c1f.fits
#iraf: 276.09, 255.68
#xx = 276.1670059634217
#yy = 255.6366143804087
#a1 = 350.2057743290616
#s1 =   2.1880983261379
#a2 = 105.0654703086852
#s2 =   4.7672184223544
#bb =  23.6716030937674
#
#x0c80106t_c1f.fits
#iraf: 292.89, 251.72
#xx = 293.0905170673360
#yy = 251.7017003991349
#a1 =  27.4652969696053
#s1 =   1.5557515660662
#a2 =  15.7672865566820
#s2 =   3.9112904777026
#bb =   4.2859960039287
#
#x0c80107t_c1f.fits
#iraf: 292.35, 250.52
#xx = 292.3438963749745
#yy = 250.6692759293081
#a1 =  31.4902885168371
#s1 =   1.6542290983903
#a2 =  11.3891475245589
#s2 =   4.5793427469043
#bb =   4.1943522516578
