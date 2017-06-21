import numpy as np
import matplotlib.pyplot as plt
from parameters import *

r = np.logspace(np.log10(3.1),np.log10(2000),600)
f = 1 - np.sqrt(3 / r)

rho1 = 9e-4 * (alpha * mbh)**(-1) * mdot**(-2) * r**(3./2) * f**(-2)
T1  = 4.9e7 * (alpha * mbh)**(-0.25) * r**(-3.0/8)
H1 = 5.5e4 * mbh * mdot * f + r*0

rho2 = 8 * (alpha*mbh)**(-0.7) * mdot**(0.4) * r**(-33./20) * f**(0.4)
T2  = 2.2e8 * (alpha*mbh)**(-0.2) * mdot**(0.4) * r**(-0.9) * f**(0.4)
H2 = 2.7e3 * alpha**(-0.1) * mbh**(0.9) * mdot**(0.2) * r**1.05 * f**0.2

rho3 = 47 * (alpha*mbh)**(-0.7) * mdot**(11.0/20) * r**(-15.0/8) * f**(11.0/20)
T3  = 6.9e7 * (alpha*mbh)**(-0.2) * mdot**(0.3) * r**(-0.75) * f**(0.3)
H3 = 1.5e3 * alpha**(-0.1) * mbh**0.9 * mdot**0.15 * r**1.125 * f**0.15

rhog = np.where(r < r12, rho1, np.where(r < r23, rho2, rho3))
Tg = np.where(r < r12, T1, np.where(r < r23, T2, T3))
Hg = np.where(r < r12, H1, np.where(r < r23, H2, H3))

fig, axes = plt.subplots(1, 3, dpi = 144, figsize = (12,4))

colors = [ '#D4532B', '#47C930', '#3670B5' ]

for a in axes:
    a.set_xscale('log')
    a.set_yscale('log')
    a.axvline(r12, color = '#AEBEB7')
    a.axvline(r23, color = '#BFD2CA')

ax = axes[0]
ax.set_title('$\\rho$')
ax.plot(r, rho1, color = colors[0])
ax.plot(r, rho2, color = colors[1])
ax.plot(r, rho3, color = colors[2])
ax.plot(r, rhog, color = '#0C0A1A')

ax = axes[1]
ax.set_title('$T$')
ax.plot(r, T1, color = colors[0])
ax.plot(r, T2, color = colors[1])
ax.plot(r, T3, color = colors[2])
ax.plot(r, Tg, color = '#0C0A1A')

ax = axes[2]
ax.set_title('$H$')
ax.plot(r, H1, color = colors[0])
ax.plot(r, H2, color = colors[1])
ax.plot(r, H3, color = colors[2])
ax.plot(r, Hg, color = '#0C0A1A')

plt.show()
