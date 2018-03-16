from par import *
from diskvert import col2python

#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#


figs = [
    ('M', '$\\dot m$', mdots_1),
    ('N', '$\\nu$', nus_1),
    ('A', '$\\alpha_{{\\rm B}}$', alphas_1),
]

#------------------------------------------------------------------------------#

cm = plt.get_cmap('coolwarm')

#------------------------------------------------------------------------------#

for sn, sl, sv in figs:

    fig, axes = plt.subplots(1, 2, figsize = (figw_lg, figw_lg / 2))

    DP = [ (i, x, col2python('data.1/01-' + ix2fn(sn,(i,)) + '.dat')) for i,x in enumerate(sv) ]

    for ax in axes.ravel():
        ax.set_xlim(0, 120)
        ax.set_yscale('log')

    for i, x, dp in DP:
        d,p = dp
        ax = axes[0]
        ax.plot(d['h'], d['temp'], color = cm(i / 8.0), label = '{} = {:.2g}'.format(sl, x))
        ax = axes[1]
        ax.plot(d['h'], d['rho'], color = cm(i / 8.0), label = '{} = {:.2g}'.format(sl, x))
    
    for ax in axes.ravel(): ax.legend(fontsize = 9)

    # ax = axes[1,0]
    # ax.set_yscale('log')
    # for i, x, dp in DP:
    #     d,p = dp
    #     ax.plot(d['z'], d['temp'])

    # plt.savefig('multi-{}.{}'.format(sn, figext))
#------------------------------------------------------------------------------#

plt.show()

#------------------------------------------------------------------------------#
