#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import numpy as np
from diskvert import *

dx,pdx = col2python('MDFX.dat')
qk,pqk = col2python('MQFK.dat')
wk,pwk = col2python('MWFK.dat')
ck,pck = col2python('MCFK.dat')
wx,pwx = col2python('MWFX.dat')
# qx,pqx = col2python('MQFX.dat')
cx,pcx = col2python('MCFX.dat')

#------------------------------------------------------------------------------#

def PL(ax, y = lambda d,p: d['temp']):
    ax.plot(dx['h'], y(dx,pdx), color = '#B0BFC2', label = 'dyfu')

    ax.plot(ck['h'], y(ck,pck), color = '#004E9D',
        linewidth = 1.5, label = 'exact RK4')
    ax.plot(wk['h'], y(wk,pwk), color = '#7DD8FF', label = 'compt RK4')
    ax.plot(qk['h'], y(qk,pqk), color = '#0079F7', label = 'compt mod RK4')

    ax.plot(wx['h'], y(wx,pwx), linestyle = '--', color = '#FF959A',
        label = 'compt X')
    ax.plot(cx['h'], y(cx,pcx), linestyle = '--', color = '#330002',
        linewidth = 1.5, label = 'exact X')

    ax.legend(loc = 'best', fontsize = 9)

#------------------------------------------------------------------------------#

fig,axes = plt.subplots(1, 2, figsize = (12,6), dpi = 118)

for ax in axes: ax.set_xlim(6,16)

axes[0].set_yscale('log')
axes[0].set_ylim(1.5e6,3.4e7)
PL(axes[0], lambda d,p: d['temp'])

def u(d,p):
    t = d['trad'] / d['temp']
    return 3 - (t + t**2 + t**3)

PL(axes[1], u)
axes[1].set_ylim(0,3)

plt.savefig('porown.png', dpi = 144)
plt.show()
