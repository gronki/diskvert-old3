# coding: utf-8
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from diskvert import *

n = 48
fnames = [ '{}.{:03d}.tar.gz'.format('MWFK',i+1) for i in range(n) ]
fnames2 = [ '{}.{:03d}.tar.gz'.format('MCFX',i+1) for i in range(n) ]


for i,fn,fn2 in zip(range(len(fnames)),fnames,fnames2):
    fig = plt.figure()

    plt.xscale('log')
    plt.xlim(1e4,1e-2)
    plt.axvline(2.0/3.0, color='#CBCCCA')

    plt.yscale('log')
    plt.ylim(1e-9,1)
    plt.ylim(1e6,1e8)

    d,p = col2python(fn)
    d2,p2 = col2python(fn2)
    #plt.title('$r$ = {:.2f}'.format(p.radius))
    plt.plot(d['tau'], d['temp'], color = '#555555', linewidth = 1.0)
    plt.plot(d2['tau'], d2['temp'], color = 'black', linewidth = 1.25)
    # plt.plot(d['tau'], d['temp'], label='$\\beta$ = {:.2f}'.format(p.beta_0), \
    #     color = cm.inferno(float(i) / float(len(fnames)-1)))

    plt.savefig('{}.T{:03d}.png'.format('F',i+1))
    plt.close(fig)

# plt.legend()
# plt.show()
