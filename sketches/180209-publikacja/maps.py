import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogLocator
from symresults import *
from itertools import product as cartproduct

#------------------------------------------------------------------------------#

from diskvert import *
from par import *

#------------------------------------------------------------------------------#

contw = lambda n: np.logspace(np.log10(1.0), np.log10(1.65), n)
LogLvl = lambda z, n: np.logspace(np.log10(np.min(z)), np.log10(np.max(z)), n)

#------------------------------------------------------------------------------#

for ds, yl, yy, ys, xl, xx, xs in dsets:

    print ds

    dset3 = read2Ddataset(lambda i,j: 'data.0/' + ix2fn(ds, (i,j)) + '.txt',
        len(yy), len(xx))

    #--------------------------------------------------------------------------#

    fig, axes = plt.subplots(2, 4, figsize = (figw_lg, figw_lg * 0.52),
        sharex = True, sharey = True)

    for ax in axes.ravel():
        ax.set_xscale(xs)
        ax.set_yscale(ys)
        ax.set_xlim(np.min(xx), np.max(xx))
        ax.set_ylim(np.min(yy), np.max(yy))

    for ax in axes[-1,:]: ax.set_xlabel(xl)
    for ax in axes[:,0]:  ax.set_ylabel(yl)

    #--------------------------------------------------------------------------#

    ax = axes[0,0]
    ax.set_title('$T_{{\\rm avg}}$ [keV]')
    cs = ax.contourf(xx, yy, DPA(dset3, 'tavg_tmin') / 11.6e6,
        np.logspace(-0.5, 1, 13), norm = LogNorm(),
        cmap = 'CMRmap') # extend = 'both'
    plt.colorbar(cs, ax = ax, ticks = [0.1, 0.3, 1, 3, 10])

    cs2 = ax.contour(xx, yy, DPA(dset3, 'taues_tmin'), 11,
        linewidths = contw(11), colors = 'white')
    ax.clabel(cs2, inline = True, fontsize = 8.5, color = 'white', fmt = '%.1f')

    #--------------------------------------------------------------------------#

    # ax = axes[1,0]
    # ax.set_title('$T_{{\\rm avg}}$ [keV] (hot)')
    # cs = ax.contourf(xx, yy, DPA(dset3, 'tavg_hard') / 11.6e6,
    #     np.logspace(-0.5,1.5,16), cmap = 'CMRmap', norm = LogNorm())
    # plt.colorbar(cs, ax = ax, ticks = [1,3,10,30,100], format = '%.1f')

    #--------------------------------------------------------------------------#

    ax = axes[0,3]
    ax.set_title('$\\chi = F_{{\\rm rad}}^{{\\rm cor}} / F_{{\\rm rad}}^{{\\rm tot}}$ and $x$')
    cs = ax.contourf(xx, yy, DPA(dset3, 'chi_tmin'),
        levels = np.linspace(0, 1, 19), cmap = 'Spectral_r')
    plt.colorbar(cs, ax = ax, ticks = [0, 0.25, 0.5, 0.75, 1])

    lvl = [1.0, 1.25, 1.50, 1.75, 2.0, 2.25, 2.50, 3.0, 3.5, 4.0, 5.0]
    cs = ax.contour(xx, yy, DPA(dset3, 'xcor'), lvl,
        linewidths = contw(len(lvl)), colors = 'black')
    ax.clabel(cs, inline = True, fontsize = 8.5, color = 'black', fmt = '%.3g')

    #--------------------------------------------------------------------------#

    ax = axes[0,1]
    ax.set_title('$\\beta$ at $\\tau = \\tau_{{\\rm cor}}$ and $\\beta_0$')
    cs = ax.contourf(xx, yy, DPA(dset3, 'beta_cor'),
         cmap = 'RdBu_r', norm = LogNorm(), levels = np.logspace(-4.0,-1.0,17))
    plt.colorbar(cs, ax = ax, ticks = [1e-4, 1e-3, 1e-2, 1e-1])

    # lvl = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 5.0, 10.0, 20.0, 35.0, 50.0, 100.0, 200.0]
    # levels = lvl,    linewidths = (np.log10(lvl) + 3.0) / 6.0,
    beta_0 = DPA(dset3, 'beta_0')
    lvl = LogLvl(beta_0, 11)
    cs2 = ax.contour(xx, yy, beta_0, lvl, linewidths = contw(len(lvl)), colors = 'black')
    ax.clabel(cs2, inline = True, fontsize = 8.5, color = 'black', fmt = '$\\beta_0$ = %.2g')

    #--------------------------------------------------------------------------#

    ax = axes[1,0]

    ax.set_title('$Y_{{\\rm avg}}$ (at $\\tau* = 1$)')
    cs = ax.contourf(xx, yy, DPA(dset3, 'compy_therm'), 9,
        cmap = 'magma', vmin = 0)
    plt.colorbar(cs, ax = ax)

    # cs2 = ax.contour(xx, yy, DPA(dset3, 'tavg_therm') / 11.6e6, 15, colors = 'white')
    # ax.clabel(cs2, inline = True, fontsize = 8, color = 'white', fmt = '%.2g')

    #--------------------------------------------------------------------------#

    ax = axes[0,2]

    ax.set_title(u'$\\rho$ at $\\tau = 1$ and $\\log \\Sigma$')

    cs = ax.contourf(xx, yy, np.log10(DPA(dset3, 'rho_phot_nh')),
        np.linspace(16, 20, 4 * 3 + 1), cmap = 'afmhot_r')
    plt.colorbar(cs, ax = ax, ticks = [16, 17, 18, 19, 20])

    cs2 = ax.contour(xx, yy, np.log10(DPA(dset3, 'coldens') / cgs_mhydr), 11, linewidths = contw(11), colors = 'black')
    ax.clabel(cs2, inline = True, fontsize = 8, fmt = '%.3g')

    #--------------------------------------------------------------------------#

    ax = axes[1,3]

    ax.set_title('$F_{{\\rm rad}} / F_{{\\rm acc}}$ and $H_{{\\rm disk}}$')

    cs = ax.contourf(xx, yy, DPA(dset3, 'fbfrac_top'), 9, cmap = 'PuBuGn', vmax = 1)
    plt.colorbar(cs, ax = ax, ticks = np.linspace(0, 1, 20 + 1))

    zz = DPA(dset3, 'hdisk')
    cs2 = ax.contour(xx, yy, zz, np.logspace(0, np.log10(np.max(zz)), 13), linewidths = contw(13), colors = 'white')
    ax.clabel(cs2, inline = True, fontsize = 8, fmt = '%.2g')

    #--------------------------------------------------------------------------#

    ax = axes[1,1]

    ax.set_title('$\\epsilon$ and $\\tau_{{\\rm es}}$ at $\\tau* = 1$')
    kes = DPA(dset3, 'kapsct_therm')
    kabp = DPA(dset3, 'kapabp_therm')

    zz = kabp / (kes + kabp)
    cs = ax.contourf(xx, yy, zz, np.logspace(-4.0,-0.5,13), norm = LogNorm(),  cmap = 'RdYlGn')
    plt.colorbar(cs, ax = ax, ticks = [1e-4, 1e-3, 1e-2, 1e-1])

    zz = DPA(dset3, 'taues_therm')
    cs = ax.contour(xx, yy, zz, np.logspace(np.log10(np.min(zz)), np.log10(np.max(zz)), 11),
        linewidths = contw(11), colors = 'black')
    ax.clabel(cs, inline = True, fontsize = 8.5, color = 'black', fmt = '%.2g')

    #--------------------------------------------------------------------------#

    ax = axes[1,2]
    ax.set_title('$\\min \\left\\{ d \\ln \\Lambda \\ / \\ d \\ln T \\right\\}$')
    cm = plt.cm.get_cmap('PiYG')
    cmx = lambda lvl: [ cm(x) for x in np.linspace(0, 1, len(lvl)) ]
    cs = ax.contourf(xx, yy, DPA(dset3, 'instabil'),
        levels = [-2.0, -1.0, -0.5, -0.2, -0.1, -0.05, -0.02, -0.01, -0.005, 0, 0.001, 0.002, 0.005, 0.010, 0.020, 0.050],
        colors = [cm(x) for x in [0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.47, 0.53, 0.6, 0.67, 0.8, 0.9, 0.95, 1.0]],
        extend = 'both')
    cs2 = ax.contour(xx, yy, DPA(dset3, 'instabil'), levels = [0], linewidths = 0.8, colors = 'black', alpha = 0.4, linestyles = ':')
    # ax.clabel(cs, fontsize = 9)
    cb = plt.colorbar(cs, ax = ax, ticks = [-1, -0.2, -0.05, -0.01, 0, 0.002, 0.01, 0.05])
    cb.add_lines(cs2)
    # cs = ax.contour(xx, yy, DPA(dset3, 'taues_therm'), 11,
    #     linewidths = np.linspace(0.2, 1.8, 11),
    #     colors = 'black')
    # ax.clabel(cs, inline = True, fontsize = 9, color = 'black', fmt = '%.1f')

    #--------------------------------------------------------------------------#

    plt.tight_layout()
    plt.subplots_adjust(wspace = 0.12, hspace = 0.17)
    plt.savefig('maps-{}.{}'.format(ds, figext))
    plt.close(fig)
