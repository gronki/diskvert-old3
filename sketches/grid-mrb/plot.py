# coding: utf-8
from diskvert import *
import matplotlib.pyplot as plt
from sys import argv
for fn in argv[1:]:
    dd,pd = col2python('{}D.dat'.format(fn))
    dw,pw = col2python('{}W.dat'.format(fn))
    dc,pc = col2python('{}C.dat'.format(fn))

    fig,axes = plt.subplots(2,3, figsize=(16,10))
    
    # ustawiamy osi X
    
    for ax in axes[0,:]:
        ax.set_xscale('linear')
        ax.set_xlim(0,60)
        ax.set_xlabel('z / H')
    for ax in axes[1,:]:
        ax.set_xscale('log')
        ax.set_xlim(1e5,1e-3)
        ax.set_xlabel('taues')
        ax.grid()
    
    # rysowanie jednej kolumny
    
    def PL(axes, y = lambda d: d['temp'], label = 'temperature'):
        
        ax = axes[0]
        ax.plot(dd['h'], y(dd), color = '#cccccc', label = u'DYF')
        ax.plot(dw['h'], y(dw), color = '#2277EE', label = u'COMPT')
        ax.plot(dc['h'], y(dc), color = 'black', label = u'COR')
        
        ax = axes[1]
        ax.plot(dd['tau'], y(dd), color = '#cccccc', label = u'DYF')
        ax.plot(dw['tau'], y(dw), color = '#2277EE', label = u'COMPT')
        ax.plot(dc['tau'], y(dc), color = 'black', label = u'COR')
        
        for ax in axes: ax.set_ylabel(label)
        

    # --- temperatura ---
    
    PL(axes[:,0], lambda d: d['temp'], u'temperature')
    
    for ax in axes[:,0]: 
        ax.set_ylim(1e5,1e9)
        ax.set_yscale('log')
        ax.legend(loc = 'upper left', fontsize = 9)
    
    # --- gestosc ---
        
    PL(axes[:,1], lambda d: d['rho'], u'density')
    
    for ax in axes[:,1]: 
        ax.set_ylim(1e-8,1)
        ax.set_yscale('log')
        ax.legend(loc = 'upper right', fontsize = 9)
    
    
    # --- beta ---
        
    PL(axes[:,2], lambda d: d['pgas'] / d['pmag'], u'plasma beta')
    
    for ax in axes[:,2]: 
        ax.set_ylim(1e4,1e-4)
        ax.set_yscale('log')
        ax.legend(loc = 'lower right', fontsize = 9)
    
    fig.suptitle('mdot = {:.3e}, r = {:.3e}, beta0 = {:.3e}'.format(pd.mdot,pd.radius,pd.beta_0))
    
    plt.subplots_adjust(0.08, 0.08, 0.96, 0.92, 0.21, 0.25)
    plt.savefig('{}.png'.format(fn), dpi = 144)
    #plt.show()
    plt.close(fig)
