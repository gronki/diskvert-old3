import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc as mplrc
from parameters import *
from numpy.linalg import solve
from sys import argv

model = 'new' if len(argv) == 1 else argv[1]

assert model in ['old','new']

calc_r123 = calc_r123_old if model == 'old' else calc_r123_new
calc_rTH = calc_rTH_old if model == 'old' else calc_rTH_new

r12,r23 = calc_r123(mbh,mdot,alpha)
r = np.logspace(np.log10(3.05),np.log10(r23**1.5),2048)
(rho1, rho2, rho3), (T1, T2, T3), (H1, H2, H3) \
    = calc_rTH(mbh,mdot,r,alpha)


rhog = np.where(r < r12, rho1, np.where(r < r23, rho2, rho3))
Tg = np.where(r < r12, T1, np.where(r < r23, T2, T3))
Hg = np.where(r < r12, H1, np.where(r < r23, H2, H3))
rhoi = np.array(rhog)
Ti = np.array(Tg)
Hi = np.array(Hg)

def ramp_smooth(it, niter, ramp0 = 0):
    t = np.clip((it + 1.0) / niter, 0, 1)
    return ramp0 + (1 - ramp0) * (3 * t**2 - 2 * t**3)

from coefficients import compcoeffs

for it in range(16):
    rmp = ramp_smooth(it,16,0.05)
    for i in range(len(r)):
        A = np.ndarray(3, order = 'F')
        M = np.ndarray((3,3), order = 'F')
        A,M = compcoeffs(mbh, mdot, r[i], alpha, rhoi[i], Ti[i], Hi[i], A, M)
        dY = solve(M,-A)
        rhoi[i] += dY[0] * rmp
        Ti[i] += dY[1] * rmp
        Hi[i] += dY[2] * rmp

A11,A12,A13 = fun(mbh,mdot,r,alpha,rho1,T1,H1)
A21,A22,A23 = fun(mbh,mdot,r,alpha,rho2,T2,H2)
A31,A32,A33 = fun(mbh,mdot,r,alpha,rho3,T3,H3)
Ag1,Ag2,Ag3 = fun(mbh,mdot,r,alpha,rhog,Tg,Hg)
Ai1,Ai2,Ai3 = fun(mbh,mdot,r,alpha,rhoi,Ti,Hi)

mplrc('font', family = 'serif', size = 10)

fig, axes = plt.subplots(2, 4, figsize = (18, 9), dpi = 118)
fig.subplots_adjust(
    left    = 0.06,
    bottom  = 0.09,
    right   = 0.97,
    top     = 0.92,
    wspace  = 0.25,
    hspace  = 0.25,
)

colors = [ '#D56323', '#47C930', '#3670B5' ]
color_g = '#0C0A1A'
color_i = '#BB23AB'

def axp(a, ys = 'log'):
    a.set_xscale('log')
    a.set_xlabel('$r / r_{\\rm schw}$')
    a.set_yscale(ys)
    if (r12 >= 3): a.axvline(r12, color = '#AEBEB7')
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

f = lambda r: np.where(r > 3, 1 - np.sqrt(3/r), 0)
lgs = lambda x: np.sign(x) * np.log10(1 + np.abs(x))

ax = axp(axes[1,0], 'linear')
ax.set_title('$A_1$ (accretion rate)')
A0 = 1e-4 * mdot * mbh * sol_mdot_edd
plot5(ax, lgs(A11/A0), lgs(A21/A0), lgs(A31/A0), lgs(Ag1/A0), lgs(Ai1/A0))

ax = axp(axes[1,1], 'linear')
ax.set_title('$A_2$ (accretion power)')
# A0 = 1e-10 * mbh**2 * mdot * cgs_facc_sun
A0 = 1e20 * alpha**(-0.8)*mbh**(-0.8)*mdot**1.6
plot5(ax, lgs(A12/A0), lgs(A22/A0), lgs(A32/A0), lgs(Ag2/A0), lgs(Ai2/A0))

ax = axp(axes[1,2], 'linear')
ax.set_title('$A_3$ (gas pressure)')
A0 = 1e12 * (alpha*mbh)**(-0.9) * mdot**(0.8)
plot5(ax, lgs(A13/A0), lgs(A23/A0), lgs(A33/A0), lgs(Ag3/A0), lgs(Ai3/A0))

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

plt.savefig('porown_{}.png'.format(model))
plt.show()
