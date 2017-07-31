# coding: utf-8
import matplotlib.pyplot as plt
import numpy as np
from sys import argv
from diskvert import col2python

d0,p0 = col2python('MCFK.dat')
d1,p1 = col2python('MDFX.dat')

name = argv[1]

frames = [ int(x) for x in argv[2:] ]

for f in frames:
    fig, axes = plt.subplots(1, 2, figsize = (9.0,4.5))

    axes[0].plot(d1['h'], d1['temp'], linewidth = 1.0, color = '#37B916', label = 'dyfu')
    axes[0].plot(d0['h'], d0['temp'], linewidth = 1.5, color = '#02060E', label = 'RK4')
    axes[1].plot(d1['h'], d1['rho'], linewidth = 1.0, color = '#37B916', label = 'dyfu')
    axes[1].plot(d0['h'], d0['rho'], linewidth = 1.5, color = '#02060E', label = 'RK4')

    axes[0].set_yscale('log')
    axes[0].set_ylim(1e6,1e9)
    axes[1].set_yscale('log')
    axes[1].set_ylim(1e-8,1)

    d,p = col2python('{n}.dat[{f}]'.format(f = f, n = name))
    axes[0].plot(d['h'], d['temp'], label = 'relax', color='#363636')
    axes[1].plot(d['h'], d['rho'], label = 'relax', color='#363636')

    axes[0].legend(fontsize = 9)
    axes[1].legend(fontsize = 9)

    plt.savefig('{n}.T{f:03d}.png'.format(f = f, n = name), dpi = 144)
    plt.close(fig)
