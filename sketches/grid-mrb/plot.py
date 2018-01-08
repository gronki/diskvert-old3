# coding: utf-8
from diskvert import *
import matplotlib.pyplot as plt
from sys import argv

for fn in argv[1:]:
    dd,pd = col2python('data/{}D.dat'.format(fn))
    dw,pw = col2python('data/{}W.dat'.format(fn))
    dc,pc = col2python('data/{}C.dat'.format(fn))

    fig,axes = plt.subplots(2,3, figsize=(14,9))

    clr = dict(
        d = '#AAAAAA',
        w = '#134DE0',
        c = '#171512',
    )

    # --------------

    ax = axes[0,0]
    ax.set_title(u'temperature')

    ax.set_xlim(0,60)

    ax.set_yscale('log')
    ax.set_ylim(1e5,1e9)

    ax.plot(dd['h'], dd['temp'], color = clr['d'], label = 'DYF')
    ax.plot(dw['h'], dw['temp'], color = clr['w'], label = 'COMPT')
    ax.plot(dc['h'], dc['temp'], color = clr['c'], label = 'BIL')
    ax.legend(fontsize = 9, loc = 'best')
    ax.grid()

    # --------------

    ax = axes[1,0]
    ax.set_title(u'temperature')

    ax.set_xlim(1e4,1e-4)
    ax.set_xscale('log')

    ax.set_yscale('log')
    ax.set_ylim(1e5,1e9)

    ax.plot(dd['taues'], dd['temp'], color = clr['d'], label = 'DYF')
    ax.plot(dw['taues'], dw['temp'], color = clr['w'], label = 'COMPT')
    ax.plot(dc['taues'], dc['temp'], color = clr['c'], label = 'BIL')
    ax.grid()

    ax.axvline(pw.taues_cor, linestyle = '--', color = clr['w'])
    ax.axvline(pc.taues_cor, linestyle = '--', color = clr['c'])
    ax.axhline(pw.tcor, linewidth = 0.5, color = clr['w'])
    ax.axhline(pc.tcor, linewidth = 0.5, color = clr['c'])

    # --------------

    ax = axes[0,1]
    ax.set_title(u'pressure')

    ax.set_xlim(0,60)

    ax.set_yscale('log')
    ax.set_ylim(1e7,1e18)

    ax.plot(dd['h'], dd['pgas'], color = clr['d'])
    ax.plot(dw['h'], dw['pgas'], color = clr['w'])
    ax.plot(dc['h'], dc['pgas'], color = clr['c'], label = 'gas')
    ax.plot(dd['h'], dd['prad'], ':', color = clr['d'])
    ax.plot(dw['h'], dw['prad'], ':', color = clr['w'])
    ax.plot(dc['h'], dc['prad'], ':', color = clr['c'], label = 'rad')
    ax.plot(dd['h'], dd['pmag'], '--', color = clr['d'])
    ax.plot(dw['h'], dw['pmag'], '--', color = clr['w'])
    ax.plot(dc['h'], dc['pmag'], '--', color = clr['c'], label = 'mag')
    ax.legend(fontsize = 9, loc = 'best')
    ax.grid()

    # --------------

    ax = axes[1,1]
    ax.set_title(u'beta')
    ax.set_xlim(0,60)
    ax.set_yscale('log')
    ax.set_ylim(1e3,1e-4)

    ax.plot(dd['h'], dd['beta'], color = clr['d'])
    ax.plot(dw['h'], dw['beta'], color = clr['w'])
    ax.plot(dc['h'], dc['beta'], color = clr['c'], label = 'mag')
    ax.plot(dd['h'], dd['pgas'] / dd['prad'], '--', color = clr['d'])
    ax.plot(dw['h'], dw['pgas'] / dw['prad'], '--', color = clr['w'])
    ax.plot(dc['h'], dc['pgas'] / dc['prad'], '--', color = clr['c'], label = 'rad')
    ax.legend(fontsize = 9, loc = 'lower right')
    ax.grid()

    # --------------

    ax = axes[0,2]
    ax.set_title(u'flux')
    ax.set_xlim(0,60)
    ax.set_ylim(0,1)

    ax.plot(dd['h'], dd['frad'] / pd.facc, color = clr['d'])
    ax.plot(dw['h'], dw['frad'] / pd.facc, color = clr['w'])
    ax.plot(dc['h'], dc['frad'] / pd.facc, color = clr['c'], label = 'rad')
    ax.plot(dd['h'], dd['fmag'] / pd.facc, '--', color = clr['d'])
    ax.plot(dw['h'], dw['fmag'] / pd.facc, '--', color = clr['w'])
    ax.plot(dc['h'], dc['fmag'] / pd.facc, '--', color = clr['c'], label = 'mag')
    ax.legend(fontsize = 9, loc = 'center right')
    ax.grid()

    # --------------

    ax = axes[1,2]
    ax.set_title(u'beta')

    # ax.set_xlim(0,60)
    ax.set_xlim(1e4,1e-4)
    ax.set_xscale('log')

    ax.set_yscale('log')
    ax.set_ylim(1e3,1e-4)

    ax.plot(dd['taues'], dd['beta'], color = clr['d'])
    ax.plot(dw['taues'], dw['beta'], color = clr['w'])
    ax.plot(dc['taues'], dc['beta'], color = clr['c'], label = 'mag')
    ax.plot(dd['taues'], dd['pgas'] / dd['prad'], '--', color = clr['d'])
    ax.plot(dw['taues'], dw['pgas'] / dw['prad'], '--', color = clr['w'])
    ax.plot(dc['taues'], dc['pgas'] / dc['prad'], '--', color = clr['c'], label = 'rad')
    ax.legend(fontsize = 9, loc = 'lower right')
    ax.grid()

    # --------------

    from numpy import log10
    fig.suptitle('mdot = {:.4f}, r = {:.2f}, a = {:.5f}, beta0 = {}'.format(pd.mdot, pd.radius, pd.alpha, pd.beta_0))

    plt.subplots_adjust(0.08, 0.08, 0.96, 0.92, 0.21, 0.25)
    plt.savefig('img/{}.png'.format(fn))
    plt.close(fig)
