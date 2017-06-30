# coding: utf-8
from diskvert import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

mbh, mdot, alpha = 10, 0.01, 0.1

r = np.logspace(0.5, 3, 2**8)
d = np.logspace(-2,1, 2**7)

rho = np.ndarray((d.shape[0], r.shape[0]), order = 'F')
T = np.ndarray((d.shape[0], r.shape[0]), order = 'F')

apxdisk2d(mbh, mdot, r, alpha, d, rho, T, r.shape[0], d.shape[0])

fig,axes = plt.subplots(2, 1, figsize = (12,8))
plt.subplots_adjust(0.09, 0.09, 0.98, 0.92, 0.20, 0.27)
ax = axes[0]
m = ax.pcolor(r, d, rho, cmap = 'hot')
# ax.set_aspect(1)
ax.loglog()
ax.set_title('Density [${\\rm g} / {\\rm cm}^3$]')
plt.colorbar(m, ax = ax)

ax = axes[1]
m = ax.pcolor(r, d, T, cmap = 'hot')
# ax.set_aspect(1)
ax.loglog()
ax.set_title('Temperature [$\\rm K$]')
plt.colorbar(m, ax = ax)

# plt.savefig('dysk2d.png')
plt.show()
