# coding: utf-8

import matplotlib.pyplot as plt

from diskvert import col2python
d0,p0 = col2python('cond.dat')

for i in range(p0.niter):
    d,p = col2python('cond.dat[{}]'.format(i))
    fig, axes = plt.subplots(1,2,figsize=(9,4))
    ax = axes[0]
    ax.set_yscale('log')
    ax.set_ylim(1e6,1e10)
    ax.plot(d['h'], d['trad'], color = '#A2A2A2')
    ax.plot(d['h'], d['temp'], color = '#222222')

    ax = axes[1]
    ax.plot(d['h'], d['frad'], color = '#FF8A49')
    ax.plot(d['h'], d['fmag'], color = '#2E84F6')
    ax.plot(d['h'], d['fcnd'], color = '#14EE2A')
    ax.plot(d['h'], d['frad'] + d['fmag'] + d['fcnd'], color = '#000000', linestyle = '--')
    plt.savefig("img.{:03d}.png".format(i))
