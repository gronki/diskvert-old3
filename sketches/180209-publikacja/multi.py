from par import *
from diskvert import col2python
from diskvert.cgs import *

#------------------------------------------------------------------------------#

cm = plt.get_cmap('coolwarm')

#------------------------------------------------------------------------------#

for sn, sl, sv in dsets_multi:

    fig, axes = plt.subplots(2, 3, figsize = (figw_lg, 0.5 * figw_lg), sharex = True)

    DP = [ (i, x, col2python('data.1/multi-' + ix2fn(sn,(i,)) + '.dat')) for i,x in enumerate(sv) ]

    for ax in axes.ravel():
        ax.set_xlim(1e4,1e-4)
        ax.axvline(1, linestyle = ':', color = '#B6BDAB')
        ax.set_xscale('log')

    for ax in axes[-1,:]:
        ax.set_xlabel('$\\tau$')

    for i, x, dp in DP:
        d,p = dp
        xx = d['tau'] #/ p.zdisk
        clr = cm(i / float(len(sv) - 1))
        lbl = '{} = {:.2g}'.format(sl, x)

        ax = axes[0,0]
        ax.plot(xx, d['temp'], color = clr, label = lbl)
        ax.set_yscale('log')
        ax.set_title('$T_{{\\rm gas}}$')

        ax = axes[0,1]
        ax.plot(xx, d['rho'] / cgs_mhydr, color = clr, label = lbl)
        ax.set_yscale('log')
        ax.set_title('$\\rho$')

        ax = axes[1,0]
        ax.plot(xx, d['beta'], color = clr, label = lbl)
        # ax.plot(xx, (d['pgas'] + d['prad']) / d['pmag'], color = clr, label = lbl)
        ax.set_ylim(1e1,1e-4)
        ax.set_yscale('log')
        ax.set_title('$\\beta$')

        ax = axes[1,1]
        ax.plot(xx, d['frad'] / p.facc, color = clr, label = lbl)
        ax.set_title('$F_{{\\rm rad}} / F_{{\\rm acc}}$')

        ax = axes[0,2]
        ax.plot(xx, d['pmag'], color = clr, label = lbl)
        ax.set_yscale('log')
        ax.set_title('$P_{{\\rm mag}}$')

        ax = axes[1,2]
        ax.plot(xx, d['kabp'] / d['ksct'], color = clr, label = lbl)
        ax.set_yscale('log')
        ax.set_title('$\\epsilon$')


    li, la = axes[0,0].get_legend_handles_labels()
    fig.legend(li, la, fontsize = 9, loc = 'center right')

    plt.tight_layout()
    plt.subplots_adjust(right = 0.87)
    plt.savefig('multi-{}.{}'.format(sn, figext))
#------------------------------------------------------------------------------#

# plt.show()

#------------------------------------------------------------------------------#
