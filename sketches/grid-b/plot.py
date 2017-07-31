# coding: utf-8
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from diskvert import *

fnames = [ 'MCFK.{0:04d}.tar.gz'.format(i+1) for i in range(16) ]


for i,fn in zip(range(len(fnames)),fnames):
    fig = plt.figure()

    plt.xscale('log')
    plt.xlim(1e4,1e-2)
    plt.axvline(2.0/3.0, color='#CBCCCA')

    plt.yscale('log')
    plt.ylim(1e-9,1)
    plt.ylim(1e6,1e8)

    d,p = col2python(fn)
    plt.title('$\\beta$ = {:.2f}'.format(p.beta_0))
    plt.plot(d['tau'], d['temp'], color = 'black', linewidth = 1.25)
    # plt.plot(d['tau'], d['temp'], label='$\\beta$ = {:.2f}'.format(p.beta_0), \
    #     color = cm.inferno(float(i) / float(len(fnames)-1)))

    plt.savefig('T.{0:04d}.png'.format(i+1))
    plt.close(fig)

# plt.legend()
# plt.show()
