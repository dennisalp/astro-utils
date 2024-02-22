import sys

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sts

dat= np.loadtxt(sys.argv[1])
xx=np.linspace(0,100,1000001)
aa,bb=np.histogram(dat[:,2], xx)
nn=np.sum(aa)
cc=1-np.cumsum(aa)/nn

da2= np.loadtxt(sys.argv[2])
a2,b2=np.histogram(da2[:,2], xx)
n2=np.sum(a2)
c2=1-np.cumsum(a2)/n2

plt.semilogy(xx[:-1] , cc, lw=3)
plt.semilogy(xx[:-1] , c2)
plt.semilogy(xx[:-1], sts.chi2.sf(xx[:-1], 2))
plt.show()
