from pdb import set_trace as st
from glob import glob
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits



def open_fits(ff):
    dat = fits.open(ff)
    dat = dat[0].data
    dat = np.nan_to_num(dat)
    return dat



of = open_fits
bb = sorted(glob('*f438w*'))
rr = sorted(glob('*f625w*'))

bb = of(bb[1]) - (of(bb[0]) + of(bb[2])) / 2
rr = of(rr[1]) - (of(rr[0]) + of(rr[2])) / 2

bb = bb[565:610, 570:615]
rr = rr[565:610, 570:615]


plt.imshow(bb, origin='lower', vmin=np.percentile(bb,10), vmax=np.percentile(bb,90))
plt.title('f438w')
plt.figure()

plt.imshow(rr, origin='lower', vmin=np.percentile(rr,10), vmax=np.percentile(rr,90))
plt.title('f625w')
plt.show()

st()
