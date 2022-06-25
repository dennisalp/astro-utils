import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import subprocess
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

ndata = 65
xfit, yfit = 2.10451, 0.57855
xmin, xmax = 1.95, 2.25
ymin, ymax = 0.45, 0.7
fit_chi = 264.82
x_label = '$\\Gamma$'
y_label = '$f_\\mathrm{pc}$'
out_name = 'fc_gamma'

x = np.linspace(xmin, xmax, ndata)
y = np.linspace(ymin, ymax, ndata)
Z = np.loadtxt(out_name+'.dat')
Z = Z-fit_chi
levels = np.array([2.3,4.61,9.21]) # Delta fit statistic, 68, 90 and 99%
X,Y=np.meshgrid(x,y)

plt.figure(figsize=(6, 4.5), dpi=300)
CS = plt.contour(X, Y, Z, levels=levels, colors=('r', 'g', 'b'))

fmt = {}
strs = ['68~\%','90~\%','99~\%']
for l, s in zip(CS.levels, strs):
    fmt[l] = s
plt.clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=13)

im = plt.imshow(Z, interpolation='bilinear', origin='lower', extent=(xmin, xmax, ymin, ymax), norm=colors.LogNorm(vmin=Z.min(), vmax=Z.max()), cmap=cm.gray, aspect='auto')
# make a colorbar for the contour lines
CB = plt.colorbar(im, extend='max')
CB.ax.set_ylabel('$\\Delta\\chi^2$')

plt.plot(xfit, yfit, '.', color='orange')

plt.xlabel(x_label)
plt.ylabel(y_label)

plt.savefig(out_name+'.pdf',bbox_inches='tight', pad_inches=0.03)
#subprocess.call(['pdfcrop',out_name+'.pdf'], shell=False)
#subprocess.call(['rm',out_name+'.pdf'], shell=False)
##subprocess.call(['mkdir','-p','plots/'+path_stem], shell=False)
#subprocess.call(['mv',out_name+'-crop.pdf',out_name+'.pdf'], shell=False)
##plt.show()
