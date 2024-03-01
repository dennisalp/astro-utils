import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import VBBinaryLensing

mpl.rc('font',**{'family':'serif','serif':['Palatino']})
mpl.rc('text', usetex=True)
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 20

# Initialize VBBinaryLensing() class object, set relative accuracy
VBBL = VBBinaryLensing.VBBinaryLensing()
VBBL.RelTol = 1e-03
t  =  np.linspace(7100, 7200, 500)

s = 1.2 # separation between the two lenses in units of total ang. Einstein radii
q = 0.1 # mass ratio: mass of the lens on the right divided by mass of the lens on the left
rho = 0.01 # source radius in Einstein radii of the total mass.
alpha = 0. # angle between lens axis and source trajectory
tE = 100.3 # einstein radius crossing time
t0 = 7154. # time of peak magnification
u0 = 0.3 # impact parameter

# Position of the center of the source with respect to the center of mass.
tau = (t - t0)/tE
y1 = -u0*np.sin(alpha) + tau*np.cos(alpha)
y2 = u0*np.cos(alpha) + tau*np.sin(alpha)

mag = np.zeros(len(tau))


params = [np.log(s), np.log(q), u0, alpha, np.log(rho), np.log(tE), t0]
results = VBBL.BinaryLightCurve(params, t)
mag = results[0]
y1 = results[1] #Note that source coordinates are calculated by BinaryLightCurve itself and are given as output
y2 = results[2]

# Calculate the cirtical curves and the caustic curves
caustics = VBBL.Caustics(s, q)    # Returns caustics only (much simpler use in Python)

fig, ax = plt.subplots(figsize=(12,7))

ax.plot(t, mag, 'k-')
ax.grid(True)

ax2 = fig.add_axes([.54, .44, .4, .4], aspect=1)
for caustic in caustics:
    ax2.plot(caustic[0], caustic[1], 'k-')
ax2.plot(y1, y2, 'k--')
ax2.grid(True)
ax2.set_xlim(-1,1)
ax2.set_ylim(-1,1)
plt.show()
