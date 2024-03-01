from pdb import set_trace as st

import numpy as np
from scipy.optimize import minimize, brute
import matplotlib.pyplot as plt
import scipy.stats as sts

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

import VBBinaryLensing



################################################################
# Functions
def plt_lcs():
    plt.figure()
    plt.plot(bb, cc/dt, 'k.')
    plt.plot(b1[:-1], c1/np.diff(b1)[0], drawstyle='steps-post')
    plt.plot(b2[:-1], c2/np.diff(b2)[0], drawstyle='steps-post')
    plt.plot(b3[:-1], c3/np.diff(b3)[0], drawstyle='steps-post')

    plt.grid()
    plt.xlabel('$t$ (s)')
    plt.ylabel('Count rate (s$^{-1})$')
    plt.show()


def get_mag(params):
    results = mod.BinaryLightCurve(params, tt)
    mag = np.array(results[0])
    y1 = np.array(results[1])
    y2 = np.array(results[2])
    return mag, y1, y2


def get_logl(xx):

    amp = xx[0]
    params = xx[1:]

    mag, _, _ = get_mag(params)
    logl = -np.sum(sts.poisson(mag*amp).logpmf(c1))
    print_pars(xx, logl)
    return logl


def print_pars(xx, logl=None):
    if logl and logl > 10000:
        return
    elif logl:
        buf = 'log(L): {0:9.2f}, '.format(logl)
    else:
        buf = ''
    buf += 'amp: {0:8.4f}, '
    buf += 'ss: {1:8.4f}, '
    buf += 'qq: {2:8.4f}, '
    buf += 'u0: {3:8.4f}, '
    buf += 'aa: {4:8.4f}, '
    buf += 'rr: {5:8.4f}, '
    buf += 'te: {6:8.4f}, '
    buf += 't0: {7:8.4f}'
    print(buf.format(
        xx[0],
        np.exp(xx[1]),
        np.exp(xx[2]),
        xx[3],
        xx[4],
        np.exp(xx[5]),
        np.exp(xx[6]),
        xx[7],
    ))




# generator function iterating over _sols, _curve, or _point objects 
# making use of the next keyword
def iterate_from(item):
    while item is not None:
        yield item
        item = item.next


def plt_solutions(tt, mag, caustics, criticals, amp, y1, y2):

    fig, ax = plt.subplots(figsize=(12,7))

    ax.plot(bb, cc/dt, 'k.')
    ax.plot(b1[:-1], c1/np.diff(b1)[0], color='gray', drawstyle='steps-post')

    ax.plot(tt*spd, amp*mag, 'k-')
    ax.set_xlabel('$t$ (s)')
    ax.set_ylabel('$A$')
    ax.grid(True)

    fig, ax2 = plt.subplots(figsize=(12,7))

    for ii, caustic in enumerate(caustics):
        ax2.plot(caustic[0], caustic[1], '-', color='tab:blue', label=f'Caustic {ii}')
    for ii, critical in enumerate(criticals):
        ax2.plot(critical[0], critical[1], '-', color='tab:orange', label=f'Critical Curve {ii}')

    ax2.plot(y1, y2, 'k.', ms=1, label='Source path')
    ax2.legend()
    ax2.grid(True)
    # ax2.set_xlim(2*min(y1.min(), -y1.max()), 2*max(-y1.min(), y1.max()))
    # ax2.set_ylim(2*min(y2.min(), -y2.max()), 2*max(-y2.min(), y2.max()))

    plt.show()



################################################################
# Parameters
spd = 24*3600
tt = np.loadtxt('tt.csv')
tt = tt[115:2080]
t0 = np.median(tt)
tt -= t0
dt = np.diff(np.unique(tt)).min()



################################################################
# Time series
bb, cc = np.unique(tt, return_counts=True)
b1 = np.arange(tt.min(), tt.max()+dt, dt)
c1, _ = np.histogram(tt, b1)
b2 = np.arange(tt.min(), tt.max()+dt, dt*10)
c2, _ = np.histogram(tt, b2)
b3 = np.arange(tt.min(), tt.max()+dt, dt*100)
c3, _ = np.histogram(tt, b3)

assert np.all(np.diff(tt) >= 0), 'events must be chronological'

# plt_lcs()



################################################################
# Microlens modeling
mod = VBBinaryLensing.VBBinaryLensing()
mod.RelTol = 1e-03
# mod.Tol = 1e-2

tt = (b1[1:] + b1[:-1]) / 2
tt /= spd

if False:
    amp = 2
    ss = 1.2  # separation between the two lenses in units of total ang. Einstein radii
    qq = .1  # mass ratio: mass of the lens on the right divided by mass of the lens on the left
    u0 = 0.3  # impact parameter
    aa = 0  # angle between lens axis and source trajectory
    rr = .02  # source radius in Einstein radii of the total mass.
    te = .002  # einstein radius crossing time
    t0 = .0  # time of peak magnification
    params = [np.log(ss), np.log(qq), u0, aa, np.log(rr), np.log(te), t0]
elif True:
    amp = 0.201
    ss = 1.5  # separation between the two lenses in units of total ang. Einstein radii
    qq = 1.  # mass ratio: mass of the lens on the right divided by mass of the lens on the left
    u0 = -0.  # impact parameter
    aa = 0.  # angle between lens axis and source trajectory
    rr = .002  # source radius in Einstein radii of the total mass.
    te = .005  # einstein radius crossing time
    t0 = .0  # time of peak magnification
    params = [np.log(ss), np.log(qq), u0, aa, np.log(rr), np.log(te), t0]
elif False:
    amp = 1
    ss = 0.4  # separation between the two lenses in units of total ang. Einstein radii
    qq = 0.01  # mass ratio: mass of the lens on the right divided by mass of the lens on the left
    u0 = 1.206  # impact parameter
    aa = -2.75  # angle between lens axis and source trajectory
    rr = .0002  # source radius in Einstein radii of the total mass.
    te = .2  # einstein radius crossing time
    t0 = .3456  # time of peak magnification
    params = [np.log(ss), np.log(qq), u0, aa, np.log(rr), np.log(te), t0]
elif True:
    amp = 0.0201
    ss = 0.8298  # separation between the two lenses in units of total ang. Einstein radii
    qq = 0.0324  # mass ratio: mass of the lens on the right divided by mass of the lens on the left
    u0 = -0.1418  # impact parameter
    aa = 1.1461  # angle between lens axis and source trajectory
    rr = .00002  # source radius in Einstein radii of the total mass.
    te = .3285  # einstein radius crossing time
    t0 = .0198  # time of peak magnification

    xx = [amp, np.log(ss), np.log(qq), u0, aa, np.log(rr), np.log(te), t0]
    fit = minimize(get_logl, xx, options={'maxiter': 300})
    print_pars(fit.x, logl=fit.fun)
    amp = fit.x[0]
    params = fit.x[1:]
elif True:
    amp = slice(0.0002, 0.03, 0.004)
    ss = [np.log(0.1), np.log(5)]  # separation between the two lenses in units of total ang. Einstein radii
    qq = [np.log(.001), np.log(.02)]  # mass ratio: mass of the lens on the right divided by mass of the lens on the left
    u0 = [-.05, .05]  # impact parameter
    aa = [-3, 3]  # angle between lens axis and source trajectory
    rr = [np.log(.0001), np.log(.001)]  # source radius in Einstein radii of the total mass.
    te = [np.log(.012), np.log(0.25)]  # einstein radius crossing time
    t0 = slice(0,1,1)  # time of peak magnification

    ranges = [amp, ss, qq, u0, aa, rr, te, t0]
    fit = brute(get_logl, ranges, Ns=5)
    print_pars(fit)
    amp = fit[0]
    params = fit[1:]
else:
    amp = [0.0003, 0.03]
    ss = [np.log(0.1), np.log(5)]  # separation between the two lenses in units of total ang. Einstein radii
    qq = [np.log(.001), np.log(.02)]  # mass ratio: mass of the lens on the right divided by mass of the lens on the left
    u0 = [-.1, .1]  # impact parameter
    aa = [.3, 2]  # angle between lens axis and source trajectory
    rr = [np.log(.0001), np.log(.001)]  # source radius in Einstein radii of the total mass.
    te = [np.log(.012), np.log(0.2)]  # einstein radius crossing time
    t0 = [-.0003, 0.]  # time of peak magnification

    ranges = [amp, ss, qq, u0, aa, rr, te, t0]
    fit = brute(get_logl, ranges, Ns=5)
    print_pars(fit)
    amp = fit[0]
    params = fit[1:]


ss = np.exp(params[0])
qq = np.exp(params[1])
mag, y1, y2 = get_mag(params)
caustics = mod.Caustics(ss, qq) # Returns _sols object containing n crit. curves followed by n caustic curves
criticals = mod.CriticalCurves(ss, qq)
plt_solutions(tt, mag, caustics, criticals, amp, y1, y2)

# st()
