python -c "import sys
import numpy as np
from astropy.cosmology import FlatLambdaCDM
zz = float(sys.argv[1])

mpc = 3.086e24

cos = FlatLambdaCDM(H0=70, Om0=0.27, Tcmb0=2.725)
mu = cos.distmod(zz).value
dl = cos.luminosity_distance(zz).value
dc = cos.comoving_distance(zz).value
da = cos.angular_diameter_distance(zz).value
scale = 2*np.pi/(360*3600)*da*1e3
f2l = 4*np.pi*(dl*mpc)**2
l2f = f2l**-1

print('{0:13.6f} Mpc (luminosity)'.format(dl))
print('{0:13.6f} Mpc (comoving)'.format(dc))
print('{0:13.6f} Mpc (angular)'.format(da))
print('{0:13.6f} kpc/\"'.format(scale))
print('{0:13.2e} (flux to luminosity)'.format(f2l))
print('{0:13.2e} (luminosity to flux)'.format(l2f))
print('{0:13.6f} mag (distmod)'.format(mu))
" ${1}
