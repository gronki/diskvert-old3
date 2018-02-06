from sys import argv
import numpy as np
import matplotlib.pyplot as plt
from diskvert import *
import re

#------------------------------------------------------------------------------#

datext = r'\.(dat|tar(\.[gx]z|)|t[gx]z)$'

for fn in argv[1:]:

    #--------------------------------------------------------------------------#

    d,p = col2python(fn)
    if not p.converged: print u'WARNING: {} not converged'.format(fn)

    def setx(ax):
        ax.set_xlabel('z / H')
        ax.set_xscale('log')
        # ax.set_xlim(1, 4e2)
        ax.set_xlim(2, None)

    fig, axes = plt.subplots(2, 3, figsize = (12.5, 8.5), sharex = True)
    fig.canvas.set_window_title(fn)
    bigtitle = u'$\\log M_{{\\rm BH}}$ = {:.1f}, $\\dot{{m}}$ = {:.5f}, $R / R_{{\\rm schw}}$ = {:.1f}, $\\alpha_{{\\rm B}}$ = {:.3f}, $\\eta$ = {:.3f}, $\\nu_{{\\rm rec}}$ = {:.1f}'.format(np.log10(p.mbh), p.mdot, p.radius, p.alpha, p.eta, p.nu)
    plt.suptitle(bigtitle)

    #--------------------------------------------------------------------------#

    ax = axes[0,0]
    ax.set_title(u'Temperature')

    ax.plot(d['h'], d['trad'], label = '$T_{\\rm rad}$', color = '#B6BBBF')
    ax.plot(d['h'], d['tavg'], label = '$T_{\\rm avg}$', color = '#ABC897', linestyle = '--')
    ax.plot(d['h'], d['temp'], label = '$T_{\\rm gas}$', color = '#DB4024')

    ax.axvline(p.zphot  / p.zscale,  linewidth = 0.7, color = '#4664D0', label = '$\\tau = 1$')
    ax.axvline(p.ztherm / p.zscale, linewidth = 0.7, color = '#F14628', label = '$\\tau* = 1$')
    if hasattr(p,'zeqbc'):
        ax.axvline(p.zeqbc / p.zscale, linewidth = 0.7, color = '#55C222', label = '$\\Lambda_B / \\Lambda_C = 1$')
    if p.has_corona:
        ax.axhline(p.tavg_tmin, color = '#ABC897', linewidth = 0.7, label = 'T average')

    setx(ax)

    ax.set_yscale('log')
    ax.legend(loc = 'best', fontsize = 9)


    #--------------------------------------------------------------------------#

    ax = axes[0,1]
    ax.set_title(u'Pressure')

    ax.plot(d['h'], d['prad'], label = '$P_{\\rm rad}$', color = '#4EBD12')
    ax.plot(d['h'], d['pmag'], label = '$P_{\\rm mag}$', color = '#2C8BED')
    ax.plot(d['h'], d['pgas'], label = '$P_{\\rm gas}$', color = '#DB4024')
    ax.loglog()

    setx(ax)


    ax.legend(loc = 'best', fontsize = 9)

    #--------------------------------------------------------------------------#

    ax = axes[1,0]
    ax.set_title(u'Energy transport')

    ax.plot(d['h'], d['frad'] / p.facc, label = '$F_{\\rm rad}$', color = '#4EBD12')
    ax.plot(d['h'], d['fmag'] / p.facc, label = '$P_{\\rm mag}$', color = '#2C8BED')
    ax.plot(d['h'], d['heat'] * d['z'] / p.facc, color = '#F069D7', label = '$z \\cdot \\Gamma_{{\\rm mag}}$')
    ax.set_ylim(0, 1)

    setx(ax)

    ax.axvline(p.zphot  / p.zscale,  linewidth = 0.7, color = '#4664D0')
    ax.axvline(p.ztherm / p.zscale, linewidth = 0.7, color = '#F14628')
    if p.has_corona and hasattr(p,'zeqbc'):
        ax.axvline(p.zeqbc / p.zscale, linewidth = 0.7, color = '#55C222')


    ax.legend(loc = 'best', fontsize = 9)

    #--------------------------------------------------------------------------#

    ax = axes[1,1]
    ax.set_title(u'Heating and cooling')

    ax.plot(d['h'], d['heat'], color = '#941B84', label = 'heating $\\Gamma_{{\\rm mag}}$')
    ax.plot(d['h'], d['coolb'], '--', color = '#F49425', label = 'brehms.')
    ax.plot(d['h'], d['coolc'], ':', color = '#17CF67', label = 'compton')
    ax.set_ylim(0, max(d['heat']) * 1.2)

    ax.axvline(p.zphot  / p.zscale,  linewidth = 0.7, color = '#4664D0')
    ax.axvline(p.ztherm / p.zscale, linewidth = 0.7, color = '#F14628')
    if p.has_corona and hasattr(p,'zeqbc'):
        ax.axvline(p.zeqbc / p.zscale, linewidth = 0.7, color = '#55C222')

    ax2 = ax.twinx()
    ax2.plot(d['h'], d['compfr'], color = '#A1ABB0', label = '$\\Lambda_C / (\\Lambda_B + \\Lambda_C)$')
    ax2.set_ylim(0,1)

    # ax.set_yscale('log')
    Li1, La1 = ax.get_legend_handles_labels()
    Li2, La2 = ax2.get_legend_handles_labels()
    ax.legend(Li1 + Li2, La1 + La2, fontsize = 9, loc = 'best')

    setx(ax)


    #--------------------------------------------------------------------------#

    ax = axes[1,2]
    ax.set_title(u'Dimensionless variables')

    ax.plot(d['h'], d['taues'], color = '#AAB5B0', label = '$\\tau_{{\\rm es}}$')
    ax.plot(d['h'], d['beta'], color = '#245ECD', label = '$\\beta_{{\\rm mag}}$')
    ax.set_yscale('log')
    ax.set_ylim(1e4,1e-4)

    ax2 = ax.twinx()
    # ax2.set_ylim(1.2,1.8)
    # ax2.axhline(4.0 / 3, linestyle = '--', linewidth = 0.7, color = '#F1C0B5')
    # ax2.axhline(5.0 / 3, linestyle = '--', linewidth = 0.7, color = '#B5F1B5')
    # ax2.plot(d['h'], d['adiab1'], color = '#EC691D', label = '$\\gamma$')
    ax2.plot(d['h'], d['compy'], label = u'total Y', color = '#D314E6')
    ax2.set_ylim(0,0.25)

    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines + lines2, labels + labels2, fontsize = 9, loc = 'best')

    ax.axvline(p.zphot  / p.zscale,  linewidth = 0.7, color = '#4664D0')
    ax.axvline(p.ztherm / p.zscale, linewidth = 0.7, color = '#F14628')
    if p.has_corona and hasattr(p,'zeqbc'):
        ax.axvline(p.zeqbc / p.zscale, linewidth = 0.7, color = '#55C222')

    setx(ax)

    #--------------------------------------------------------------------------#

    ax = axes[0,2]
    ax.set_title(u'Density')
    ax.plot(d['h'], d['rho'], color = '#357100')
    ax.set_yscale('log')

    ax.axvline(p.zphot  / p.zscale,  linewidth = 0.7, color = '#4664D0')
    ax.axvline(p.ztherm / p.zscale, linewidth = 0.7, color = '#F14628')
    if p.has_corona and hasattr(p,'zeqbc'):
        ax.axvline(p.zeqbc / p.zscale, linewidth = 0.7, color = '#55C222')

    setx(ax)

    #--------------------------------------------------------------------------#

    plt.subplots_adjust(0.05, 0.08, 0.95, 0.90, 0.27, 0.21)

    if re.search(datext, fn):
        plt.savefig(re.sub(datext, '.png', fn))
    else:
        plt.savefig(fn + '.png')

    #--------------------------------------------------------------------------#

# plt.show()

#------------------------------------------------------------------------------#
