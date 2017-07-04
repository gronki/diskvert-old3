from sys import argv
import numpy as np
import matplotlib.pyplot as plt
from diskvert import *

d = np.loadtxt(argv[1])

z = d[:,1]
rho = d[:,2]
T = d[:,3]
Trad = d[:,4]
Frad = d[:,5]
Prad = 4 * cgs_stef / (3 * cgs_c) * Trad**4

fig, axes = plt.subplots(2,2)

ax = axes[0,0]
ax.plot(z,rho)

ax = axes[0,1]
ax.plot(z,Trad)
ax.plot(z,T)

kappa = lambda rho,T: rho * T**(-3.5)

ax = axes[1,0]
ax.plot(z, kappa(rho,T))
ax.set_yscale('log')

ax = axes[1,1]

ax.plot(z, ( 4 * cgs_stef * T**4 - 3 * cgs_c * Prad ) * rho * kappa(rho,T) \
        + 12 * cgs_kapes * rho * Prad * cgs_boltz * ( T - Trad ) \
        / ( cgs_mel * cgs_c ))

plt.show()
