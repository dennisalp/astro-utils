from pdb import set_trace as st
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from glob import glob



x0, y0 = 4070, 4115
bins = np.linspace(5000,10000,100)
files = {
    'acisf24295N001_evt2.fits': (),
    'acisf21304N002_evt2.fits': (),
    'acisf24654N001_evt2.fits': (),
    'acisf22849N002_evt2.fits': (),
    'acisf25514N001_evt2.fits': (),
    'acisf23534N002_evt2.fits': (),
    'acisf22425N001_evt2.fits': (),
    'acisf25906N001_evt2.fits': (),
    'acisf24652N001_evt2.fits': (),
}


# align
for zoom in np.arange(1,2.5,0.7):
    for ff, xy in files.items():
        dat = fits.open(ff)[1].data
        xx = dat['x']
        yy = dat['y']
        ee = dat['energy']

        xy = xy if xy else (x0, y0)

        ii = (np.abs(xx - xy[0]) < 15/zoom) & (np.abs(yy - xy[1]) < 20/zoom)
        xy = (xx[ii].mean(), yy[ii].mean())
        files[ff] = xy

    #     plt.figure()
    #     plt.plot(xx[ii], yy[ii], '.', ms=1)
    #     plt.plot(xy[0], xy[1], 'ok', ms=5)
    #     plt.grid()
    #     plt.title(ff)

    # plt.show()


xx = np.array([], dtype=np.float32)
yy = np.array([], dtype=np.float32)
ee = np.array([], dtype=np.float32)
gr = np.array([], dtype=np.float32)
# stack
for ff, xy in files.items():
    dat = fits.open(ff)[1].data
    xx = np.concatenate((xx, dat['x'] - xy[0]))
    yy = np.concatenate((yy, dat['y'] - xy[1]))
    ee = np.concatenate((ee, dat['energy']))
    gr = np.concatenate((gr, dat['grade']))

ii = (np.abs(xx) < 10) & (np.abs(yy) < 10)
fe = ii & (ee > 6400) & (ee <  6900) & (gr <= 6) & (gr != 1) & (gr != 5)
nt = ii & (ee > 7400) & (ee < 10000) & (gr <= 6) & (gr != 1) & (gr != 5)

plt.figure()
plt.hist(ee[ii], bins=bins, color='gray')
plt.hist(ee[fe], bins=bins, color='r')
plt.hist(ee[nt], bins=bins, color='b')
plt.grid()

plt.figure()
plt.plot(xx[fe], yy[fe], '.r', ms=1)
plt.plot(xx[nt], yy[nt], '.b', ms=1)
plt.plot(0, 0, 'ok', ms=5)
plt.grid()
plt.show()

st()
