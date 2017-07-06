import numpy as np
import matplotlib.pyplot as plt
from diskvert import *

def name2description(name):
    from re import match

    m = match(r"(M|A)(D|W|C|M)(F?B?Z?E?)(C?)(K|X)",name)
    mm = dict(
        heat = m.group(1),
        temp = m.group(2),
        opacity = m.group(3),
        condu = m.group(4) == 'C',
        solver = m.group(5),
    )

    tempmet = dict(
        D = 'dyfu',
        W = 'compton only',
        C = 'corona',
        M = 'multi temp',
    )
    heatmet = dict(
        M = 'magnetic',
        A = 'alpha',
    )
    opamet = dict(
        F = 'ff',
        B = 'bf',
        Z = 'bb',
        E = 'es',
    )
    solmet = dict(K = 'RK4', X = 'relax')

    return " ".join([
        solmet[mm['solver']],
        heatmet[mm['heat']], tempmet[mm['temp']],
        "+".join([opamet[x] for x in mm['opacity']]) \
    ] + (['condu'] if 'C' in mm['opacity'] else []))

def plot(name, z, tau, rho, temp, flux, ztop):

    longname = name2description(name)

    fig, axes = plt.subplots(1,3,figsize=(12,4))
    plt.subplots_adjust(0.07, 0.15, 0.97, 0.88, 0.24, 0.26)
    fig.suptitle(longname)

    axes[0].plot(tau, rho)
    axes[0].set_ylabel('density')
    axes[0].set_yscale('log')
    axes[0].set_ylim(3e-9,1)
    axes[1].plot(tau, temp)
    axes[1].set_ylabel('temperature')
    axes[1].set_ylim(0,2e7)
    axes[2].plot(tau, flux)
    axes[2].set_ylabel('flux rad')
    axes[2].set_ylim(0,4e20)

    for ax in axes:
        ax.set_xscale('log')
        ax.set_xlim(1e4,1e-4)
        ax.set_xlabel('tau')

    plt.savefig('T{}.png'.format(name), dpi = 144)
    plt.close(fig)

    # --------------------------------------------------

    fig, axes = plt.subplots(1,3,figsize=(12,4))
    plt.subplots_adjust(0.07, 0.15, 0.97, 0.88, 0.24, 0.26)
    fig.suptitle(longname)

    axes[0].plot(z * 1e-5, rho)
    axes[0].set_ylabel('density')
    axes[0].set_yscale('log')
    axes[0].set_ylim(3e-9,1e0)
    axes[1].plot(z * 1e-5, temp)
    axes[1].set_ylabel('temperature')
    axes[1].set_ylim(0,2e7)
    axes[2].plot(z * 1e-5, flux)
    axes[2].set_ylabel('flux rad')
    axes[2].set_ylim(0,4e20)

    for ax in axes:
        ax.set_xlim(0,ztop * 1e-5)
        ax.set_xlabel('z [km]')

    plt.savefig('Z{}.png'.format(name), dpi = 144)
    plt.close(fig)

from filelist import files

for fn,ft in files.items():
    if ft == 'AK':
        d,p = col2python('{}.dat'.format(fn))
        plot(fn, d['z'], d['tau_es'], d['rho'], d['tgas'], d['frad'], \
                12 * p.zscale)
    if ft == 'AX':
        d,p = col2python(fn+'.dat')
        plot(fn, d['z'], d['tau'], d['rho'], d['temp'], d['frad'], \
                12 * p.zscale)
    if ft == 'MK':
        d,p = col2python(fn+'.dat')
        plot(fn, d['z'], d['tau'], d['rho'], d['temp'], d['frad'], \
                80 * p.zscale)
    if ft == 'MX':
        d,p = col2python(fn+'.dat')
        plot(fn,  d['z'], d['tau'], d['rho'], d['temp'], d['frad'], \
                80 * p.zscale)
