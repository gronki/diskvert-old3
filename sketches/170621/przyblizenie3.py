import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc as mplrc
from parameters import *
from numpy.linalg import solve

def fun2(r,rho,T,H):
    A = np.ndarray(3)
    M = np.ndarray((3,3))
    A[0] = 4*H**3*alpha*cgs_pi*rho*sqrt(cgs_graw*cgs_msun/(cgs_rschw_sun \
      **3*mbh**2*r**3)) - cgs_mcrit_sun*mbh*mdot*(-1.73205080756888* \
      sqrt(1.0/r) + 1)
    A[1] = 0.5625*H**4*alpha*rho**2*(cgs_graw*cgs_msun/(cgs_rschw_sun** \
          3*mbh**2*r**3))**1.5*(T**(-3.5)*kappa_0*rho + cgs_kapes) - 4* \
          T**4*cgs_stef
    A[2] = -H**2*rho*(cgs_graw*cgs_msun/(cgs_rschw_sun**3*mbh**2*r**3))** \
          1.0 + (1.0/3.0)*T**4*cgs_a + 2*T*cgs_boltz*rho/cgs_mhydr
    M[0, 0] = 4*H**3*alpha*cgs_pi*sqrt(cgs_graw*cgs_msun/(cgs_rschw_sun**3* \
          mbh**2*r**3))
    M[1, 0] = 0.5625*H**4*T**(-3.5)*alpha*kappa_0*rho**2*(cgs_graw* \
          cgs_msun/(cgs_rschw_sun**3*mbh**2*r**3))**1.5 + 1.125*H**4* \
          alpha*rho*(cgs_graw*cgs_msun/(cgs_rschw_sun**3*mbh**2*r**3))** \
          1.5*(T**(-3.5)*kappa_0*rho + cgs_kapes)
    M[2, 0] = -H**2*(cgs_graw*cgs_msun/(cgs_rschw_sun**3*mbh**2*r**3))** \
          1.0 + 2*T*cgs_boltz/cgs_mhydr
    M[0, 1] = 0
    M[1, 1] = -1.96875*H**4*T**(-4.5)*alpha*kappa_0*rho**3*(cgs_graw* \
          cgs_msun/(cgs_rschw_sun**3*mbh**2*r**3))**1.5 - 16*T**3* \
          cgs_stef
    M[2, 1] = (4.0/3.0)*T**3*cgs_a + 2*cgs_boltz*rho/cgs_mhydr
    M[0, 2] = 12*H**2*alpha*cgs_pi*rho*sqrt(cgs_graw*cgs_msun/(cgs_rschw_sun \
          **3*mbh**2*r**3))
    M[1, 2] = 2.25*H**3*alpha*rho**2*(cgs_graw*cgs_msun/(cgs_rschw_sun**3* \
          mbh**2*r**3))**1.5*(T**(-3.5)*kappa_0*rho + cgs_kapes)
    M[2, 2] = -2*H*rho*(cgs_graw*cgs_msun/(cgs_rschw_sun**3*mbh**2*r**3))** \
          1.0
    return (A,M)


r = np.logspace(np.log10(3.2),np.log10(5*r23),300)
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
rhoi = np.array(rhog)
Ti = np.array(Tg)
Hi = np.array(Hg)

def ramp_smooth(it, niter, ramp0 = 0):
    t = np.clip(it / (niter - 1.0), 0, 1)
    return ramp0 + (1 - ramp0) * (3 * t**2 - 2 * t**3)

for it in range(16):
    rmp = ramp_smooth(it,12,0.02)
    for i in range(len(r)):
        A,M = fun2(r[i], rhoi[i], Ti[i], Hi[i])
        dY = solve(M,-A)
        rhoi[i] += dY[0] * rmp
        Ti[i] += dY[1] * rmp
        Hi[i] += dY[2] * rmp

A11,A12,A13 = fun(r,rho1,T1,H1)
A21,A22,A23 = fun(r,rho2,T2,H2)
A31,A32,A33 = fun(r,rho3,T3,H3)
Ag1,Ag2,Ag3 = fun(r,rhog,Tg,Hg)
Ai1,Ai2,Ai3 = fun(r,rhoi,Ti,Hi)

mplrc('font', family = 'serif', size = 10)

format_A4_vertical = (8.267, 11.692)
format_A4_horizontal = (11.692, 8.267)

fig, axes = plt.subplots(2, 4, figsize = format_A4_horizontal)
fig.subplots_adjust(
    left    = 0.08,
    bottom  = 0.10,
    right   = 0.95,
    top     = 0.93,
    wspace  = 0.44,
    hspace  = 0.41,
)

colors = [ '#D56323', '#47C930', '#3670B5' ]
color_g = '#0C0A1A'
color_i = '#BB23AB'

def axp(a,ys = 'log'):
    a.set_xscale('log')
    a.set_xlabel('$r / r_{\\rm schw}$')
    a.set_yscale(ys)
    a.axvline(r12, color = '#AEBEB7')
    a.axvline(r23, color = '#BFD2CA')
    return a

def plot5(ax,a1,a2,a3,ag,ai):
    ax.plot(r, a1, color = colors[0], linewidth = 0.8, label = '1')
    ax.plot(r, a2, color = colors[1], linewidth = 0.8, label = '2')
    ax.plot(r, a3, color = colors[2], linewidth = 0.8, label = '3')
    ax.plot(r, ag, color = color_g, linewidth = 1.1, label = 'sklej')
    ax.plot(r, ai, color = color_i, linewidth = 1.4, label = 'relax')

ax = axp(axes[0,0])
ax.set_title('$\\rho$')
ax.set_ylabel('$\\rho$ [${\\rm g}/{\\rm cm}^3$]')
plot5(ax, rho1, rho2, rho3, rhog, rhoi)
ax.legend(fontsize = 8)

ax = axp(axes[0,1])
ax.set_title('$T$')
ax.set_ylabel('$T$ [${\\rm K}$]')
plot5(ax, T1, T2, T3, Tg, Ti)

ax = axp(axes[0,2])
ax.set_title('$H$')
ax.set_title('$H$ [$\\rm cm$]')
plot5(ax, H1, H2, H3, Hg, Hi)

ax = axp(axes[1,0], 'linear')
ax.set_title('$A_1$')
plot5(ax, A11, A21, A31, Ag1, Ai1)

ax = axp(axes[1,1], 'linear')
ax.set_title('$A_2$')
plot5(ax, A12, A22, A32, Ag2, Ai2)

ax = axp(axes[1,2], 'linear')
ax.set_title('$A_3$')
plot5(ax, A13, A23, A33, Ag3, Ai3)

ax = axp(axes[0,3])
ax.set_title('$\\kappa$')
ax.set_xscale('log')
ax.set_yscale('log')
ax.plot(r, kappa(rhog,Tg), color = color_g)
ax.plot(r, kappa(rhoi,Ti), color = color_i)

ax = axp(axes[1,3])
ax.set_title('$P_{\\rm rad} / P_{\\rm gas}$')
ax.set_xscale('log')
ax.set_yscale('log')
betarad = lambda rho,T: (cgs_a / 3 * T**4) / (2 * cgs_boltz * rho * T / cgs_mhydr)
ax.plot(r, betarad(rhog,Tg), color = color_g)
ax.plot(r, betarad(rhoi,Ti), color = color_i)

plt.savefig('porown.pdf')
plt.show()
