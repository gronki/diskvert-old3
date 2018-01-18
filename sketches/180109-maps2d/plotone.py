from sys import argv
import numpy as np
import matplotlib.pyplot as plt
from diskvert import *

#------------------------------------------------------------------------------#

d,p = col2python(argv[1])

if not p.converged: print u'not converged'

plt.figure(figsize = (10,10))
plt.suptitle('R = {radius:.1f} x Rschw, zeta = {eta:.3f}, beta0 = {beta:.2e}'.format(eta = p.eta, beta = p.beta_0, radius = p.radius))

#------------------------------------------------------------------------------#

plt.subplot(221)

plt.plot(d['h'], d['trad'], label = '$T_{\\rm rad}$', color = '#B6BBBF')
plt.plot(d['h'], d['tavg'], label = '$T_{\\rm avg}$', color = '#ABC897', linestyle = '--')
plt.plot(d['h'], d['temp'], label = '$T_{\\rm gas}$', color = '#DB4024')

plt.axvline(p.zphot  / p.zscale,  linewidth = 0.7, color = '#4664D0', label = '$\\tau = 1$')
plt.axvline(p.ztherm / p.zscale, linewidth = 0.7, color = '#F14628', label = '$\\tau* = 1$')
plt.axvline(p.zeqbc / p.zscale, linewidth = 0.7, color = '#55C222', label = '$\\Lambda_B / \\Lambda_C = 1$')

plt.axhline(p.tcor, color = '#ABC897', linewidth = 0.7, label = 'T average')

plt.xlabel('z / H')
plt.xlim(1, 2e2)
plt.xscale('symlog', linthreshx = 1)

plt.yscale('log')
plt.legend(loc = 'best', fontsize = 10)


#------------------------------------------------------------------------------#

plt.subplot(222)
plt.plot(d['h'], d['prad'], label = '$P_{\\rm rad}$', color = '#4EBD12')
plt.plot(d['h'], d['pmag'], label = '$P_{\\rm mag}$', color = '#2C8BED')
plt.plot(d['h'], d['pgas'], label = '$P_{\\rm gas}$', color = '#DB4024')
plt.loglog()

plt.xlabel('z / H')
plt.xlim(1, 2e2)

plt.legend(loc = 'best', fontsize = 10)

#------------------------------------------------------------------------------#

plt.subplot(223)
plt.plot(d['h'], d['frad'], label = '$F_{\\rm rad}$', color = '#4EBD12')
plt.plot(d['h'], d['fmag'], label = '$P_{\\rm mag}$', color = '#2C8BED')
plt.plot(d['h'], d['heat'] * d['z'], color = '#F069D7', label = '$z \\cdot \\Gamma_{{\\rm mag}}$')
plt.ylim(0, p.facc)

plt.xlim(1, 2e2)
plt.xlabel('z / H')
plt.xscale('log')

plt.axvline(p.zphot  / p.zscale,  linewidth = 0.7, color = '#4664D0', label = '$\\tau = 1$')
plt.axvline(p.ztherm / p.zscale, linewidth = 0.7, color = '#F14628', label = '$\\tau* = 1$')
plt.axvline(p.zeqbc / p.zscale, linewidth = 0.7, color = '#55C222', label = '$\\Lambda_B / \\Lambda_C = 1$')

plt.legend(loc = 'best', fontsize = 10)

#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#

plt.subplot(224)

# plt.plot(d['taues'], d['trad'], label = '$T_{\\rm rad}$', color = '#B6BBBF')
# plt.plot(d['taues'], d['tavg'], label = '$T_{\\rm avg}$', color = '#ABC897')
# plt.plot(d['taues'], d['temp'], label = '$T_{\\rm gas}$', color = '#DB4024')
#
# # plt.axvline(p.zcor   / p.zscale,   linewidth = 0.7, color = '#C32E19')
# plt.axvline(1,  linewidth = 0.7, color = '#4664D0', label = 'photosphere')
# plt.axvline(p.taues_therm, linewidth = 0.7, color = '#62B158', label = 'therm. depth')
#
# plt.axhline(p.tcor, color = '#ABC897', linewidth = 0.7, label = 'T average')
#
# plt.xlabel('taues')
#

# L1 = (9 * d['trad']**4 - d['temp']**4) * d['kabs'] + 8 * d['ksct'] * d['trad']**4 * cgs_boltz * d['temp'] / (cgs_mel * cgs_c**2)
# L2 = d['heat'] / (2 * cgs_stef * d['rho'])


# plt.plot(d['h'], d['coolc'] / d['coolb'])
plt.plot(d['h'], d['heat'], color = '#941B84', label = 'heating $\\Gamma_{{\\rm mag}}$')
plt.plot(d['h'], d['coolb'] * d['cool0'], '--', color = '#F49425', label = 'brehms.')
plt.plot(d['h'], d['coolc'] * d['cool0'], ':', color = '#17CF67', label = 'compton')
plt.ylim(min(d['heat']), None)

plt.axvline(p.zphot  / p.zscale,  linewidth = 0.7, color = '#4664D0')
plt.axvline(p.ztherm / p.zscale, linewidth = 0.7, color = '#F14628')
plt.axvline(p.zeqbc / p.zscale, linewidth = 0.7, color = '#55C222')

plt.yscale('log')
plt.xscale('log')
plt.xlim(1, 2e2)
plt.legend(loc = 'best', fontsize = 10)

#------------------------------------------------------------------------------#

plt.show()

#------------------------------------------------------------------------------#
