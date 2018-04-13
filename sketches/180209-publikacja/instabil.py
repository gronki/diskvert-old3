import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from symresults import *
from itertools import product as cartproduct

#------------------------------------------------------------------------------#

from diskvert import *
from par import *

#------------------------------------------------------------------------------#

nn5 = len(mdots_5)
dp = [ col2python('data.1/05-' + ix2fn('M', (i,)) + '.dat') for i in range(nn5) ]
zmax = 90

uli = lambda n: np.linspace(-1,1,n)
lili = lambda lo, hi, n: hi - np.abs(np.linspace(lo - hi, hi - lo, n))
lilg = lambda lo, hi, n: np.exp(lili(np.log(lo), np.log(hi), n))
es = lambda x,a: (np.sign(x) * x**2 + x * a) / (1 + a)
es2 = lambda x,a: (1 - a) * np.sign(x) * x**2 + a * x
es3 = lambda x,a: (1 - a)              * x**3 + a * x
def georan(n):
    t = np.arange(n) - 0.5 * (n - 1)
    return np.sign(t) * 2**(np.abs(t) - 0.5 * (n - 1))

#------------------------------------------------------------------------------#

fig, axes = plt.subplots(nn5, 3, figsize = (figw_lg, figw_lg * 1.15), sharex = 'all', sharey = 'row')

for it, mdot_ in enumerate(mdots_5):

    d,p = dp[it]

    Tlo = p.teff * 0.78
    Thi = p.teff * 60

    a1 = (11.0 / 3)**0.25
    a2 = (8 * p.kappa_es * cgs_boltz)**2 / (3 * p.kabp0 * cgs_mel * cgs_c**2)**2 * (1e6)**10
    a3 = 8.0 / 3.0 * (3.0 / 11.0)**0.125 * p.kappa_es / p.kabp0 * cgs_boltz / (cgs_mel * cgs_c**2) * (1e6)**4.5

    #--------------------------------------------------------------------------#

    for ax in axes[it,:]:
        ax.set_xlim(0, zmax)
        ax.set_yscale('log')
        ax.set_ylim(Tlo, Thi)

    #--------------------------------------------------------------------------#

    nx, ny = 300, 200
    T = np.logspace(np.log10(Tlo), np.log10(Thi), ny)
    h = np.linspace(0, zmax, nx)
    T0 = np.interp(h, d['h'], d['temp'])
    rho = np.interp(h, d['h'], d['rho'])
    Trad = np.interp(h, d['h'], d['trad'])
    heat = np.interp(h, d['h'], d['heat'])

    bil = np.ndarray((ny,nx))
    LT = np.ndarray((ny,nx))
    BC = np.ndarray((ny,nx))

    for i,j in cartproduct(range(ny),range(nx)):
        rho0 = rho[j] * T0[j] / T[i]
        kabp = p.kabp0 * rho0 / T[i]**3.5
        Lff = 4 * cgs_stef * rho0 * kabp * (T[i]**4 - Trad[j]**4)
        LTff = 2 * cgs_stef * rho0 * kabp / T[i] \
            * (11 * Trad[j]**4 - 3 * T[i]**4)
        Les = 4 * cgs_stef * rho0 * p.kappa_es * Trad[j]**4 \
            * cgs_boltz * (4 * T[i] - 4 * Trad[j]) / (cgs_mel * cgs_c**2)
        LTes = 16 * cgs_stef * rho0 * p.kappa_es * Trad[j]**5 / T[i] \
            * cgs_boltz / (cgs_mel * cgs_c**2)
        bil[i,j] = 1 - (Lff + Les) / heat[j]
        LT[i,j] = (LTff + LTes) * T[i] / (Lff + Les)
        BC[i,j] = Les / (Lff + Les)
        if T[i] <= Trad[j]: LT[i,j] = np.ma.masked

    #--------------------------------------------------------------------------#

    ax = axes[it,0]

    cs1 = ax.contourf(h, T, BC, levels = np.linspace(0,1,13), cmap = 'RdBu_r', extend = 'both')

    # all balance solutions and solutions for gas and radiation temp
    ax.plot(d['h'], d['trad'], color = '#EEE5C2', linewidth = 1)
    ax.contour(h, T, bil, [0], linewidths = 2.5, colors = '#071803')
    cs = ax.contour(h, T, bil, 0.2 * georan(7), linewidths = lilg(0.2, 1.5, 7), colors = '#071803')
    # ax.clabel(cs, fontsize = 8)
    ax.plot(d['h'], d['temp'], ':', color = '#AFAFAF', linewidth = 2.5)

    # ax.plot(d['h'], a1 * d['trad'], color = '#9516BB', linewidth = 1)
    # ax.plot(d['h'], a2 * (d['trad'] / 1e6)**10 / d['rho']**2, '--', color = '#9516BB', linewidth = 1)

    #--------------------------------------------------------------------------#

    ax = axes[it,1]

    cs2 = ax.contourf(h, T, bil, 0.2 * es3(uli(19), 0.3), extend = 'both', cmap = 'bil1')

    # all balance solutions and solutions for gas and radiation temp
    ax.plot(d['h'], d['trad'], color = '#B2BEBF', linewidth = 1)
    ax.contour(h, T, bil, [0], linewidths = 2.5, colors = '#071803')
    ax.plot(d['h'], d['temp'], ':', color = '#AFAFAF', linewidth = 2.5)

    # the contours of dL / dT
    ax.contour(h, T, LT, georan(7) * 0.1, linewidths = lilg(0.2, 1.3, 7), colors = '#E3FFE8')

    #--------------------------------------------------------------------------#

    ax = axes[it,2]

    cs3 = ax.contourf(h, T, LT, es3(uli(15), 0.5) * 0.2, cmap = 'PiYG', extend = 'both', intervals = 'equal')
    ax.contour(h, T, LT, [0], colors = '#070E1F', linewidths = 1.0)

    # all balance solutions and solutions for gas and radiation temp
    ax.plot(d['h'], d['trad'], color = '#B2BEBF', linewidth = 1)
    ax.contour(h, T, bil, [0], linewidths = 2.5, colors = '#071803')
    ax.plot(d['h'], d['temp'], ':', color = '#AFAFAF', linewidth = 2.5)

    # ax.plot(d['h'], a2 * (d['trad'] / 1e6)**10 / d['rho']**2, '--', color = '#CAB9CF', linewidth = 0.7)
    ax.plot(d['h'], a1 * d['trad'], '-.', color = '#F18023', linewidth = 1.3)
    ax.plot(d['h'], d['rho']**2 * d['temp']**2 / (a2 * (d['trad'] / 1e6)**10), '--', color = '#F18023', linewidth = 1.3)
    # ax.plot(d['h'], d['rho'] * d['temp'] / (a3 * (d['trad'] / 1e6)**4.5))
#------------------------------------------------------------------------------#

axes[0,0].set_title(u'$\\Lambda_C / (\\Lambda_B + \\Lambda_C)$')
axes[0,1].set_title(u'$1 - \\Lambda_{{\\rm net}} \\  / \\  \\Gamma_{{\\rm mag}}$')
axes[0,2].set_title(u'$d \\ln \\Lambda_{{\\rm net}} \\ / \\ d \\ln T$')

for ax in axes[-1,:]: ax.set_xlabel('$z / H$')
for ax in axes[:,0]: ax.set_ylabel('$T$ [K]')

for ax, mdot_ in zip(axes[:,0], mdots_5):
    ax.annotate('$\\dot{{m}} = {:.3f}$'.format(mdot_), (0.05, 0.87), xycoords = 'axes fraction', color = 'white', fontsize = 11)

plt.tight_layout()
plt.subplots_adjust(bottom = 0.076)
plt.draw()

def globar(cs, ax, **kwargs):
    p0 = ax.get_position().get_points().flatten()
    cax = fig.add_axes([p0[0], 0.020, p0[2] - p0[0], 0.012])
    fig.colorbar(cs, cax = cax, orientation = 'horizontal', **kwargs)

globar(cs1, axes[-1,0], ticks = [0.0, 0.25, 0.5, 0.75, 1.0])
globar(cs2, axes[-1,1], ticks = [-0.15, -0.06, -0.02, 0, 0.02, 0.06, 0.15])
globar(cs3, axes[-1,2], ticks = [-0.2, -0.1, -0.04, 0.0, 0.04, 0.1, 0.2])

plt.savefig('instabil.{}'.format(figext))
