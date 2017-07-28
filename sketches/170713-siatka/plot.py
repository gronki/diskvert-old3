# coding: utf-8
from diskvert import *
import matplotlib.pyplot as plt
from sys import argv
for fn in argv[1:]:
    dd,pd = col2python('{}D.dat'.format(fn))
    dw,pw = col2python('{}W.dat'.format(fn))
    dc,pc = col2python('{}C.dat'.format(fn))

    fig,axes = plt.subplots(2,3, figsize=(16,10))
    
    for ax in axes[0,:]:
        ax.set_xscale('linear')
        ax.set_xlim(0,40)
        ax.set_xlabel('asdsa')
        ax.set_ylabel('asdsa')
    for ax in axes[1,:]:
        ax.set_xscale('log')
        ax.set_xlim(1e5,1e-3)
        ax.set_xlabel('asdsa')
        ax.set_ylabel('asdsa')
        ax.grid()
    
    ax = axes[0,0]

    ax.plot(dd['h'], dd['temp'], color = '#cccccc', label = u'DYF')
    ax.plot(dw['h'], dw['temp'], color = '#2277EE', label = u'COMPT')
    ax.plot(dc['h'], dc['temp'], color = 'black', label = u'COR')
    
    ax = axes[1,0]

    ax.plot(dd['tau'], dd['temp'], color = '#cccccc', label = u'DYF')
    ax.plot(dw['tau'], dw['temp'], color = '#2277EE', label = u'COMPT')
    ax.plot(dc['tau'], dc['temp'], color = 'black', label = u'COR')

    
    for ax in axes[:,0]: 
        ax.set_ylim(1e5,1e9)
        ax.set_yscale('log')
        ax.legend(loc = 'upper left', fontsize = 9)
    
    fig.suptitle('mdot = {:.3e}, r = {:.3e}, beta0 = {:.3e}'.format(pd.mdot,pd.radius,pd.beta_0))
    
    plt.subplots_adjust(0.08, 0.08, 0.96, 0.92, 0.21, 0.25)
    plt.savefig('{}.png'.format(fn), dpi = 144)
    plt.show()
    plt.close(fig)
