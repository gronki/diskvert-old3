# coding: utf-8
from diskvert import *
import matplotlib.pyplot as plt
from sys import argv
for fn in argv[1:]:
    d,p = col2python(fn)
    fig = plt.figure()
    plt.plot(d['tau'], d['trad'], color = '#bbbbbb')
    plt.plot(d['tau'], d['temp'], color = 'black')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('taues')
    plt.ylabel('temperature')
    plt.ylim(1e4,1e10)
    plt.xlim(1e4,1e-2)
    plt.title('mdot = {:.3e}, r = {:.3e}, beta0 = {:.3e}'.format(p.mdot,p.radius,p.beta_0))
    plt.savefig(fn.replace('.dat','.T.png'), dpi = 144)
    plt.close(fig)
