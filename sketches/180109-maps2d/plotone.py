from sys import argv
import numpy as np
import matplotlib.pyplot as plt
from diskvert import *

#------------------------------------------------------------------------------#

d,p = col2python(argv[1])

if not p.converged: print u'not converged'

fig, axes = plt.subplots(2, 2, figsize = (10,10), sharex = True)
plt.suptitle('R = {radius:.1f} x Rschw, zeta = {eta:.3f}, beta0 = {beta:.2e}'.format(eta = p.eta, beta = p.beta_0, radius = p.radius))

#------------------------------------------------------------------------------#

ax = axes[0,0]

ax.plot(d['h'], d['trad'], label = '$T_{\\rm rad}$', color = '#B6BBBF')
ax.plot(d['h'], d['tavg'], label = '$T_{\\rm avg}$', color = '#ABC897', linestyle = '--')
ax.plot(d['h'], d['temp'], label = '$T_{\\rm gas}$', color = '#DB4024')

ax.axvline(p.zphot  / p.zscale,  linewidth = 0.7, color = '#4664D0', label = '$\\tau = 1$')
ax.axvline(p.ztherm / p.zscale, linewidth = 0.7, color = '#F14628', label = '$\\tau* = 1$')
ax.axvline(p.zeqbc / p.zscale, linewidth = 0.7, color = '#55C222', label = '$\\Lambda_B / \\Lambda_C = 1$')

ax.axhline(p.tavg_tmin, color = '#ABC897', linewidth = 0.7, label = 'T average')

ax.set_xlabel('z / H')
ax.set_xlim(1, 2e2)
ax.set_xscale('symlog', linthreshx = 1)

ax.set_yscale('log')
ax.legend(loc = 'best', fontsize = 10)


#------------------------------------------------------------------------------#

ax = axes[0,1]
ax.plot(d['h'], d['prad'], label = '$P_{\\rm rad}$', color = '#4EBD12')
ax.plot(d['h'], d['pmag'], label = '$P_{\\rm mag}$', color = '#2C8BED')
ax.plot(d['h'], d['pgas'], label = '$P_{\\rm gas}$', color = '#DB4024')
ax.loglog()

ax.set_xlabel('z / H')
ax.set_xlim(1, 2e2)

ax.legend(loc = 'best', fontsize = 10)

#------------------------------------------------------------------------------#

ax = axes[1,0]

ax.plot(d['h'], d['frad'], label = '$F_{\\rm rad}$', color = '#4EBD12')
ax.plot(d['h'], d['fmag'], label = '$P_{\\rm mag}$', color = '#2C8BED')
ax.plot(d['h'], d['heat'] * d['z'], color = '#F069D7', label = '$z \\cdot \\Gamma_{{\\rm mag}}$')
ax.set_ylim(0, p.facc)

ax.set_xlim(1, 2e2)
ax.set_xlabel('z / H')
ax.set_xscale('log')

ax.axvline(p.zphot  / p.zscale,  linewidth = 0.7, color = '#4664D0', label = '$\\tau = 1$')
ax.axvline(p.ztherm / p.zscale, linewidth = 0.7, color = '#F14628', label = '$\\tau* = 1$')
ax.axvline(p.zeqbc / p.zscale, linewidth = 0.7, color = '#55C222', label = '$\\Lambda_B / \\Lambda_C = 1$')

ax.legend(loc = 'best', fontsize = 10)

#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#

ax = axes[1,1]


# ax.plot(d['taues'], d['trad'], label = '$T_{\\rm rad}$', color = '#B6BBBF')
# ax.plot(d['taues'], d['tavg'], label = '$T_{\\rm avg}$', color = '#ABC897')
# ax.plot(d['taues'], d['temp'], label = '$T_{\\rm gas}$', color = '#DB4024')
#
# # ax.axvline(p.zcor   / p.zscale,   linewidth = 0.7, color = '#C32E19')
# ax.axvline(1,  linewidth = 0.7, color = '#4664D0', label = 'photosphere')
# ax.axvline(p.taues_therm, linewidth = 0.7, color = '#62B158', label = 'therm. depth')
#
# ax.axhline(p.tavg_tmin, color = '#ABC897', linewidth = 0.7, label = 'T average')
#
# ax.set_xlabel('taues')
#

ax.plot(d['h'], d['instabil'])
ax.set_ylim(-1,1)
#
#
# ax.plot(d['h'], d['heat'], color = '#941B84', label = 'heating $\\Gamma_{{\\rm mag}}$')
# ax.plot(d['h'], d['coolb'] * d['cool0'], '--', color = '#F49425', label = 'brehms.')
# ax.plot(d['h'], d['coolc'] * d['cool0'], ':', color = '#17CF67', label = 'compton')
# ax.set_ylim(min(d['heat']), None)
#
# ax.axvline(p.zphot  / p.zscale,  linewidth = 0.7, color = '#4664D0')
# ax.axvline(p.ztherm / p.zscale, linewidth = 0.7, color = '#F14628')
# ax.axvline(p.zeqbc / p.zscale, linewidth = 0.7, color = '#55C222')
#
# ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim(1, 2e2)
ax.legend(loc = 'best', fontsize = 10)

#------------------------------------------------------------------------------#

plt.show()

#------------------------------------------------------------------------------#
