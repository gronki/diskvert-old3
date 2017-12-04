#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from diskvert import *

#-------------------------------------------------------------------------------
imdot = 5
bil = 'W'
#-------------------------------------------------------------------------------
nradii = 48
#-------------------------------------------------------------------------------

def radial(im,ir,ib):
    d = np.zeros((6,len(ir)))
    for i in ir:
        fn = 'M{0:02d}R{1:02d}B{2:02d}{3}.txt'.format(im,i,ib,bil)
        try:
            p = pyminiconf.pyminiconf(open('data/' + fn,'r'))
            d[:,i-1] = [p.radius, p.tau_cor, p.chicor, p.betacor, p.teff / 1.16e7, p.tcor / 1.16e7]
        except:
            d[:,i-1] = np.nan
    return d

d1 = radial(imdot, range(1, nradii + 1),  3)
d2 = radial(imdot, range(1, nradii + 1),  8)
d3 = radial(imdot, range(1, nradii + 1), 12)
d4 = radial(imdot, range(1, nradii + 1), 15)
d5 = radial(imdot, range(1, nradii + 1), 17)
d6 = radial(imdot, range(1, nradii + 1), 18)

#-------------------------------------------------------------------------------

fig,axes = plt.subplots(2, 2, figsize = (14,9))

def prepx(ax):
    ax.set_xscale('log')
    ax.set_xlabel('radius')
    ax.set_xlim(3,100)

clr = [
    '#2F05D4',
    '#2EB7CF',
    '#55C200',
    '#C2AE00',
    '#C20C00',
    '#F24EAB',
]

lbl = [
    u'$\\beta_0 = 5.27$',
    u'$\\beta_0 = 5.27$',
    u'$\\beta_0 = 14.2$',
    u'$\\beta_0 = 38.7$',
    u'$\\beta_0 = 105$',
    u'$\\beta_0 = 284$',
]

#-------------------------------------------------------------------------------

ax = axes[0,0]
ax.set_ylim(50,0)
ax.set_ylabel('$\\tau_{\\rm es,cor}$')
ax.set_title('optical depth ($\\tau_{\\rm es,cor}$)')
prepx(ax)

ax.plot(d1[0,:], d1[1,:], color = clr[0], label = lbl[0])
ax.plot(d2[0,:], d2[1,:], color = clr[1], label = lbl[1])
ax.plot(d3[0,:], d3[1,:], color = clr[2], label = lbl[2])
ax.plot(d4[0,:], d4[1,:], color = clr[3], label = lbl[3])
ax.plot(d5[0,:], d5[1,:], color = clr[4], label = lbl[4])
ax.plot(d5[0,:], d6[1,:], color = clr[5], label = lbl[5])
ax.legend(loc = 'best')

#-------------------------------------------------------------------------------

ax = axes[0,1]
ax.set_ylim(0,1)
ax.set_ylabel('$\\chi$')
ax.set_title('energy dissipated in corona ($\\chi$)')
prepx(ax)

ax.plot(d1[0,:], d1[2,:], color = clr[0], label = lbl[0])
ax.plot(d2[0,:], d2[2,:], color = clr[1], label = lbl[1])
ax.plot(d3[0,:], d3[2,:], color = clr[2], label = lbl[2])
ax.plot(d4[0,:], d4[2,:], color = clr[3], label = lbl[3])
ax.plot(d5[0,:], d5[2,:], color = clr[4], label = lbl[4])
ax.plot(d5[0,:], d6[2,:], color = clr[5], label = lbl[5])
ax.legend(loc = 'best')

#-------------------------------------------------------------------------------

ax = axes[1,0]
ax.set_ylabel('$T_{\\rm cor,av}$')
ax.set_title('temperature of the corona ($T_{\\rm cor,av}$)')
ax.set_ylim(0,2)
prepx(ax)

ax.plot(d1[0,:], d1[5,:], color = clr[0], label = lbl[0])
ax.plot(d2[0,:], d2[5,:], color = clr[1], label = lbl[1])
ax.plot(d3[0,:], d3[5,:], color = clr[2], label = lbl[2])
ax.plot(d4[0,:], d4[5,:], color = clr[3], label = lbl[3])
ax.plot(d5[0,:], d5[5,:], color = clr[4], label = lbl[4])
ax.plot(d5[0,:], d6[5,:], color = clr[5], label = lbl[5])
ax.legend(loc = 'best')

#-------------------------------------------------------------------------------

ax = axes[1,1]
ax.set_yscale('log')
ax.set_ylabel('$\\beta{\\rm cor,av}$')
ax.set_title('corona beta ($\\beta{\\rm cor,av}$)')
ax.set_ylim(1e2,1e-3)
prepx(ax)

ax.plot(d1[0,:], d1[3,:], color = clr[0], label = lbl[0])
ax.plot(d2[0,:], d2[3,:], color = clr[1], label = lbl[1])
ax.plot(d3[0,:], d3[3,:], color = clr[2], label = lbl[2])
ax.plot(d4[0,:], d4[3,:], color = clr[3], label = lbl[3])
ax.plot(d5[0,:], d5[3,:], color = clr[4], label = lbl[4])
ax.plot(d5[0,:], d6[3,:], color = clr[5], label = lbl[5])
ax.legend(loc = 'best')

#-------------------------------------------------------------------------------

plt.savefig('cordepth.M{:02d}{}.png'.format(imdot,bil), dpi = 144)
# plt.show()
