#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
from matplotlib import rc
from numpy import logspace, meshgrid

n = 2**8
rho0 = logspace(-10,2,n)
T0 = logspace(3,9,n)
rho, T = meshgrid(rho0, T0)

kapes = 0.34
# based on: http://adsabs.harvard.edu/full/1976ApJ...210..440B
red = ( 1 + 2.7e11 * rho / T**2 ) * ( 1 + (T / 4.5e8)**0.86 )

rc('font', family = 'serif', size = 10)

fig,ax = plt.subplots(figsize = (6,5))

mp = ax.pcolor(rho0, T0, kapes / red, cmap = 'hot', vmin = 0, vmax = kapes)
ax.loglog()
ax.set_title(u"Electron scattering opacity")
ax.set_xlabel("$\\rho$ [${\\rm g} / {\\rm cm}^3$]")
ax.set_ylabel("$T$ [$\\rm K$]")
plt.colorbar(mp, label = '$\\kappa_{\\rm es}$')

plt.savefig('scattering.png', dpi = 144)
