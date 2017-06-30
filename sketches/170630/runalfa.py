from diskvert import *
from ctypes import c_int, c_char, c_double
from numpy import ndarray
import matplotlib.pyplot as plt

mbh, mdot, r, alpha = 10, 0.002, 10, 0.01


c_Tgas = 2 - 1
c_rho  = 3 - 1
c_Trad = 4 - 1
c_heat = 5 - 1
c_heat_max = 6 - 1
c_compW = 7 - 1
c_compY = 8 - 1

nx = 2**10
ny,na = alf_getn()

tempmethod('B')

x = ndarray(nx)
y = ndarray((ny,nx), order = 'F')
dy = ndarray((ny,nx), order = 'F')
a = ndarray((na,nx), order = 'F')

alf_init(c_double(mbh), c_double(mdot), c_double(r), c_double(alpha))
alf_run(x,c_int(nx),y,dy,a)
# exit()
Pgas = y[0,:]
Prad = y[1,:]
Frad = y[2,:]
tau = y[3,:]
Tgas = a[1,:]
rho = a[2,:]
Trad = a[3,:]
heat = a[4,:]
heatmax = a[5,:]
compW = a[6,:]
compY = a[7,:]
kabs = a[8,:]
ksct = a[9,:]

fig, axes = plt.subplots(2,2, figsize = (8,8))

ax = axes[0,0]
ax.plot(x, Pgas, label = 'Pgas')
ax.plot(x, Prad, label = 'Prad')

ax = axes[0,1]
ax.plot(x, Tgas, label = 'Tgas')
ax.plot(x, Trad, label = 'Trad')
ax.set_ylim(1e3,1e8)

ax = axes[1,0]
ax.plot(x, heat, label = 'heat')
ax.plot(x, heatmax, label = 'heatmax')

ax = axes[1,1]
ax.plot(x, compW, label = 'compW')
ax.plot(x, compY, label = 'compY')
ax.plot(x, Prad / Pgas, label = 'Brad')
ax.plot(x, kabs / (kabs + ksct), label = 'kabs / ktot')
ax.set_ylim(1e-5,1e5)

for axrow in axes:
    for ax in axrow:
        # ax.loglog()
        ax.set_yscale('log')
        # ax.set_xlim(1e4,1e-4)
        ax.legend(fontsize = 8)

plt.show()
