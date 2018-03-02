from sys import argv
import re
import numpy as np
import matplotlib.pyplot as plt
from diskvert import *

datext = r'\.(dat|tar(\.[gx]z|)|t[gx]z)$'

D = lambda x: (x[1:] - x[:-1])
M = lambda x: (x[1:] + x[:-1]) / 2

for fn in argv[1:]:

    if not re.search(datext, fn):
        print 'skipping {}'.format(fn)
        continue

    d,p = col2python(fn)

    fig,axes = plt.subplots(2, 2, figsize = (12,10))
    fig.canvas.set_window_title(fn)
    bigtitle = u'$\\log M_{{\\rm BH}}$ = {:.1f}, $\\dot{{m}}$ = {:.5f}, $R / R_{{\\rm schw}}$ = {:.1f}, $\\alpha_{{\\rm B}}$ = {:.3f}, $\\eta$ = {:.3f}, $\\nu_{{\\rm rec}}$ = {:.1f}'.format(np.log10(p.mbh), p.mdot, p.radius, p.alpha, p.eta, p.nu)
    plt.suptitle(bigtitle)

    ax = axes[0,0]
    ax.plot(d['h'], d['pmag'][0] * d['h']**(-1-p.xcor), color = '#0056C8', linestyle = '--')
    y = d['z']**2 / (2 * p.zdisk_ss73**2)
    b0 = 2 * p.eta / p.alpha - 1 + p.nu
    pmag_th = d['pmag'][0] * (1 + (1 + b0 - p.nu) / (2 * (1 + b0) - p.nu) * y)**(-b0 / (1 + b0 - p.nu))
    ax.plot(d['h'], pmag_th, color = '#719EDA')
    ax.plot(d['h'], d['pmag'], color = '#0056C8')
    ax.plot(d['h'], d['pgas'], color = '#61A20F')
    ax.plot(d['h'], d['prad'], color = '#F56854')
    ax.axvline(p.hdisk_ss73 * np.sqrt((4 + p.alpha * p.nu / p.eta) * (1e-7**(-2/(p.xcor + 2)) - 1)))
    ax.plot(d['h'], (d['prad'] + d['pgas']) * p.alpha / (2 * p.eta + p.alpha * (p.nu - 1)), color = '#DDDDDD')
    ax.set_xscale('log')
    ax.set_xlim(2, 400)
    ax.set_yscale('log')
    ax.set_ylim(d['pmag'][0] * 1e-6, d['pmag'][0] * 30)

    ax = axes[0,1]
    ax.plot(M(d['h']), -D(d['pmag']) / (D(d['z']) * M(d['z']) * p.omega**2), color = '#4F8FE4', linewidth = 0.9)
    ax.plot(M(d['h']), -D(pmag_th) / (D(d['z']) * M(d['z']) * p.omega**2), color = '#B1D3FF', linewidth = 0.9)
    ax.plot(M(d['h']), -D(d['prad']) / (D(d['z']) * M(d['z']) * p.omega**2), color = '#F67765', linewidth = 0.9)
    ax.plot(M(d['h']), -D(d['pgas']) / (D(d['z']) * M(d['z']) * p.omega**2), color = '#A2E053', linewidth = 0.9)

    ax.axvline(p.hdisk_ss73 * np.sqrt((4 + p.alpha * p.nu / p.eta) * (1e-7**(-2/(p.xcor + 2)) - 1)))
    ax.plot(d['h'], d['rho'][0] * d['h']**(-3-p.xcor), color = '#7FD05E', linestyle = '--')
    ax.plot(d['h'], d['rho'], color = '#133704', linewidth = 1.8)
    ax.set_xscale('log')
    ax.set_xlim(2, 400)
    ax.set_yscale('log')
    ax.set_ylim(d['rho'][0] * 1e-14, d['rho'][0] * 3)

    ax = axes[1,0]
    ax.axvline(p.hdisk_ss73 * np.sqrt((4 + p.alpha * p.nu / p.eta) * (1e-7**(-2/(p.xcor + 2)) - 1)))
    ax.plot(d['h'], d['trad'], color = '#A7A7A7', linewidth = 0.9)
    ax.plot(d['h'], d['temp'], color = '#C23400', linewidth = 1.8)
    ax.set_xscale('log')
    ax.set_xlim(2, 400)
    ax.set_yscale('log')
    ax.set_ylim(p.teff * 0.5, p.teff * 1e4)

    ax = axes[1,1]
    ax.axvline(p.hdisk_ss73 * np.sqrt((4 + p.alpha * p.nu / p.eta) * (1e-7**(-2/(p.xcor + 2)) - 1)))
    ax.plot(d['h'], d['heat'] * d['z'], color = '#F267F1', linewidth = 1.0)
    ax.plot(d['h'], d['fmag'], color = '#31AAF7', linewidth = 1.2)
    ax.plot(d['h'], d['frad'], color = '#F44D4D', linewidth = 1.4)
    ax.set_xscale('log')
    ax.set_xlim(2, 400)
    ax.set_ylim(0, p.facc)


    plt.savefig(re.sub(datext, '.2.png', fn))
