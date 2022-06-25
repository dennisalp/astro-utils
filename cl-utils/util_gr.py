import sys
import os
import pdb

import matplotlib.pyplot as plt
import numpy as np

# Parameters, local, not at infinity
MM = float(sys.argv[1]) # M_Sun
RR = float(sys.argv[2]) # km
TT = 2e6 # K
DD = 5.12 # 10 kpc
sig= 5.67e-8 # Stefan-Boltzmann constant in SI
Lsun = 3.826e33 # erg s-1
Msun = 1.989e33 # g
uu = 1.660539040e-24 # g
mti = 43.9596901*uu # g
cc = 29979245800. # cm s-1
GG = 6.67259e-8 # cm3 g-1 s-2
hh = 6.6260755e-27  # erg s
kB = 1.380658e-16 # erg K-1
k2kev = 1./11604525
kpc = 3.086e21 # cm
kev2erg = 1.60218e-9

# Equations from becker09, page 183-
Rs = 2*GG/cc**2*Msun/1e5*MM
gr = np.sqrt(1-Rs/RR)
xx = 2*GG*MM*Msun/(RR*1e5*cc**2)

# For XSPEC blackbody
LL = 4*np.pi*sig*(1e3*RR)**2*TT**4
LL = LL*1e7 # SI2cgs
LL = gr**2*LL # 2infinity
KK = LL/1e39/DD**2

tie= 4*np.pi*sig*(1e3*RR/gr)**2*1/k2kev**4*1e7/1e39/DD**2

# For XSPEC carbatm
K2 = 1/DD**2

# carbatm limit, 0.3-10 keV
car_lim1 = 10**-12.9307*4*np.pi*(DD*10*kpc)**2
# 2-10 keV
car_lim2 = 10**-13.69*4*np.pi*(DD*10*kpc)**2

########
# Compute luminosity of power law from XSPEC
def pow2lum(alpha, kk, e0, e1):
    if alpha == 2:
        res = kk*np.log(e1/e0)*kev2erg
        return res, res*4*np.pi*(DD*10*kpc)**2
    else:
        hlp = 2-alpha
        res = kk*(e1**hlp-e0**hlp)/hlp*kev2erg
        return res, res*4*np.pi*(DD*10*kpc)**2

def tem2lum(TT):
    LL = 4*np.pi*sig*(1e3*RR)**2*TT**4
    LL = LL*1e7 # SI2cgs
    LL = gr**2*LL # 2infinity
    return LL

def lum2tem(LL):
    LL = LL/gr**2
    LL = LL/1e7
    TT = (LL/(4*np.pi*sig*(1e3*RR)**2))**0.25
    return TT

# Gravitational to baryonic mass using Eq. (36) of lattimer01
def mg2mb(Mg):
    Eb = Mg*Msun*0.3*xx/(1-0.25*xx)
    return Mg+Eb/Msun

# Angstrom to keV, for soft X-ray spectra
def aa2kev(aa):
    return cc/(aa*1e-8)*hh/kev2erg

# Basic conversions
def lam2nu(lam):
    return cc/lam

def nu2lam(nu):
    return cc/nu

# Spin-down of pulsars, see my note on 87a_tim (Eqs. C.6 and C.7) and manchester05
def pdot2lum(Pd, PP):
    II = 2./5*(MM*Msun)*(RR*1e5)**2
    return 4*np.pi**2*II*Pd/PP**3

def lum2pdot(LL, PP):
    II = 2./5*(MM*Msun)*(RR*1e5)**2
    return LL*PP**3/(4*np.pi**2*II)

def plt_ppdot(Lmin, Lmax, tmin, tmax):
    xx = np.logspace(-3, 2, 100)
    II = 2./5*(MM*Msun)*(RR*1e5)**2
    plt.loglog(xx, Lmin*Lsun*xx**3/(4*np.pi**2*II))
    plt.loglog(xx, Lmax*Lsun*xx**3/(4*np.pi**2*II))
    plt.loglog(xx, xx/(2*tmin))
    plt.loglog(xx, xx/(2*tmax))
    plt.show()
