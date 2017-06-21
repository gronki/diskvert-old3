import numpy as np
import matplotlib.pyplot as plt
from parameters import *

r = np.logspace(np.log10(3.5),np.log10(3.1*r23),200)
f = 1 - np.sqrt(3 / r)

rho1 = 9e-4 * (alpha * mbh)**(-1) * mdot**(-2) * r**(3./2) * f**(-2)
T1  = 4.9e7 * (alpha * mbh)**(-0.25) * r**(-3.0/8)
H1 = 5.5e4 * mbh * mdot * f + r*0

rho2 = 8 * (alpha*mbh)**(-0.7) * mdot**(0.4) * r**(-33./20) * f**(0.4)
T2  = 2.2e8 * (alpha*mbh)**(-0.2) * mdot**(0.4) * r**(-0.9) * f**(0.4)
H2 = 2.7e3 * alpha**(-0.1) * mbh**(0.9) * mdot**(0.2) * r**1.05 * f**0.2

rho3 = 47 * (alpha*mbh)**(-0.7) * mdot**0.55 * r**(-1.875) * f**0.55
T3  = 6.9e7 * (alpha*mbh)**(-0.2) * mdot**(0.3) * r**(-0.75) * f**(0.3)
H3 = 1.5e3 * alpha**(-0.1) * mbh**0.9 * mdot**0.15 * r**1.125 * f**0.15


rhog = np.where(r < r12, rho1, np.where(r < r23, rho2, rho3))
Tg = np.where(r < r12, T1, np.where(r < r23, T2, T3))
Hg = np.where(r < r12, H1, np.where(r < r23, H2, H3))

A11,A12,A13 = fun(r,rho1,T1,H1)
A21,A22,A23 = fun(r,rho2,T2,H2)
A31,A32,A33 = fun(r,rho3,T3,H3)
Ag1,Ag2,Ag3 = fun(r,rhog,Tg,Hg)

fig, axes = plt.subplots(2, 3, dpi = 118, figsize = (12,8))

colors = [ '#D4532B', '#47C930', '#3670B5' ]

def axp(a,ys = 'log'):
    a.set_xscale('log')
    a.set_yscale(ys)
    a.axvline(r12, color = '#AEBEB7')
    a.axvline(r23, color = '#BFD2CA')
    return a

ax = axp(axes[0,0])
ax.set_title('$\\rho$')
ax.plot(r, rho1, color = colors[0])
ax.plot(r, rho2, color = colors[1])
ax.plot(r, rho3, color = colors[2])
ax.plot(r, rhog, color = '#0C0A1A')

ax = axp(axes[0,1])
ax.set_title('$T$')
ax.plot(r, T1, color = colors[0])
ax.plot(r, T2, color = colors[1])
ax.plot(r, T3, color = colors[2])
ax.plot(r, Tg, color = '#0C0A1A')

ax = axp(axes[0,2])
ax.set_title('$H$')
ax.plot(r, H1, color = colors[0])
ax.plot(r, H2, color = colors[1])
ax.plot(r, H3, color = colors[2])
ax.plot(r, Hg, color = '#0C0A1A')

ax = axp(axes[1,0], 'linear')
ax.set_title('$A_1$')
ax.plot(r, A11, color = colors[0])
ax.plot(r, A21, color = colors[1])
ax.plot(r, A31, color = colors[2])
ax.plot(r, Ag1, color = '#0C0A1A')

ax = axp(axes[1,1], 'linear')
ax.set_title('$A_2$')
ax.plot(r, A12, color = colors[0])
ax.plot(r, A22, color = colors[1])
ax.plot(r, A32, color = colors[2])
ax.plot(r, Ag2, color = '#0C0A1A')

ax = axp(axes[1,2], 'linear')
ax.set_title('$A_3$')
ax.plot(r, A13, color = colors[0])
ax.plot(r, A23, color = colors[1])
ax.plot(r, A33, color = colors[2])
ax.plot(r, Ag3, color = '#0C0A1A')

plt.show()
