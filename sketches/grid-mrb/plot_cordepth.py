# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from diskvert import *

nradii = 8
ibet = 3

def radial(im,ir,ib):
    d = np.zeros((3,len(ir)))
    for i in ir:
        fn = 'M{0:03d}R{1:03d}B{2:03d}W.txt'.format(im,i,ib)
        p = pyminiconf.pyminiconf(open(fn,'r'))
        d[:,i-1] = [p.radius, p.taucor, p.chicor]
    return d

d1 = radial(1,range(1,1+nradii),ibet)
d2 = radial(2,range(1,1+nradii),ibet)
d3 = radial(3,range(1,1+nradii),ibet)

fig,axes = plt.subplots(2,1, figsize = (6,12))

for ax in axes:
    ax.set_xscale('log')
    ax.set_xlabel('radius')

ax = axes[0]
ax.set_ylim(1,100)
ax.set_yscale('log')
ax.set_title('tau of corona ($\\tau_{\\rm es,cor}$)')
ax.set_ylabel('$\\tau_{\\rm es,cor}$')

ax.plot(d1[0,:], d1[1,:], color = '#8E1DDD')
ax.plot(d2[0,:], d2[1,:], color = '#149884')
ax.plot(d3[0,:], d3[1,:], color = '#D2480D')
ax.legend(['mdot = 1e-3','mdot = 1e-2','mdot = 1e-1'], loc = 'best')

ax = axes[1]
ax.set_ylim(0.01,1)
ax.set_yscale('log')
ax.set_ylabel('$\\chi$')
ax.set_title('energy dissipated in corona ($\\chi$)')

ax.plot(d1[0,:], d1[2,:], color = '#8E1DDD')
ax.plot(d2[0,:], d2[2,:], color = '#149884')
ax.plot(d3[0,:], d3[2,:], color = '#D2480D')
ax.legend(['mdot = 1e-3','mdot = 1e-2','mdot = 1e-1'], loc = 'best')

plt.savefig('cordepth.png', dpi = 144)
plt.show()
