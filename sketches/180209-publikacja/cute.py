from par import *
from diskvert import col2python, cgs_mhydr

d,p = col2python('data.2/cute.dat')
if not p.converged: print u'WARNING: not converged'

def setx(ax):
    ax.set_xlabel('$z / z_{{\\rm disk}}$')
    # ax.set_xscale('log')
    # ax.set_xlim(1, 4e2)
    ax.set_xlim(0, 24)

fig, axes = plt.subplots(2, 3, figsize = (figw_lg, figw_lg * 0.66), sharex = True)
fig.canvas.set_window_title('cute')
bigtitle = u'$\\log M_{{\\rm BH}}$ = {:.1f}, $\\dot{{m}}$ = {:.2g}, $R / R_{{\\rm schw}}$ = {:.2g}, $\\alpha_{{\\rm B}}$ = {:.2g}, $\\eta$ = {:.2g}, $\\zeta$ = {:.2g}'.format(np.log10(p.mbh), p.mdot, p.radius, p.alpha, p.eta, p.alpha * p.nu / 2)
plt.suptitle(bigtitle)

#--------------------------------------------------------------------------#

ax = axes[0,0]
ax.set_title(u'Temperature')

ax.plot(d['z'] / p.zdisk, d['trad'], label = '$T_{\\rm rad}$', color = '#B6BBBF')
ax.plot(d['z'] / p.zdisk, d['tavg'] / d['tau'], label = '$T_{\\rm avg}$', color = '#ABC897', linestyle = '--')
ax.plot(d['z'] / p.zdisk, d['temp'], label = '$T_{\\rm gas}$', color = '#DB4024')

ax.axvline(p.zphot  / p.zdisk,  linewidth = 0.9, color = '#4664D0', label = '$\\tau = 1$')
ax.axvline(p.ztherm / p.zdisk, linewidth = 0.9, color = '#F14628', label = '$\\tau* = 1$')

setx(ax)

ax.set_yscale('log')
ax.set_ylim(1e6, 1e9)
ax.legend(loc = 'best', fontsize = 9)


#--------------------------------------------------------------------------#

ax = axes[0,1]
ax.set_title(u'Pressure')

ax.plot(d['z'] / p.zdisk, d['prad'], label = '$P_{\\rm rad}$', color = '#4EBD12')
ax.plot(d['z'] / p.zdisk, d['pmag'], label = '$P_{\\rm mag}$', color = '#2C8BED')
ax.plot(d['z'] / p.zdisk, d['pgas'], label = '$P_{\\rm gas}$', color = '#DB4024')

ptot = d['pgas'] + d['prad'] + d['pmag']
ax.plot(d['z'] / p.zdisk, ptot, label = '$P_{\\rm tot}$', color = '#1C1C1C')

ax.set_ylim(3e-5 * max(ptot), 1.5 * max(ptot))
ax.set_yscale('log')
setx(ax)

ax.legend(loc = 'best', fontsize = 9)

#--------------------------------------------------------------------------#

ax = axes[1,0]
ax.set_title(u'Energy transport')

ax.plot(d['z'] / p.zdisk, d['frad'] / p.facc, label = '$F_{\\rm rad}$', color = '#4EBD12')
ax.plot(d['z'] / p.zdisk, d['fmag'] / p.facc, label = '$P_{\\rm mag}$', color = '#2C8BED')
ax.plot(d['z'] / p.zdisk, d['heat'] * d['z'] / p.facc, color = '#F069D7', label = '$z \\cdot {{\\cal H}}$')
ax.set_ylim(0, 1)

setx(ax)

ax.axvline(p.zphot  / p.zdisk,  linewidth = 0.9, color = '#4664D0')
ax.axvline(p.ztherm / p.zdisk, linewidth = 0.9, color = '#F14628')

ax.legend(loc = 'best', fontsize = 9)

#--------------------------------------------------------------------------#

ax = axes[1,1]
ax.set_title(u'Heating and cooling')

ax.plot(d['z'] / p.zdisk, d['coolb'], color = '#F49425', label = 'brehms.', linewidth = 1.5)
ax.plot(d['z'] / p.zdisk, d['coolc'], color = '#17CF67', label = 'compton', linewidth = 1.5)
ax.plot(d['z'] / p.zdisk, d['heat'],  color = '#212121', label = 'heating $\\Gamma_{{\\rm mag}}$', linewidth = 1.8)
ax.set_yscale('log')
ax.set_ylim(max(d['heat']) * 1e-3, max(d['heat']) * 1.5)

ax.axvline(p.zphot  / p.zdisk,  linewidth = 0.9, color = '#4664D0')
ax.axvline(p.ztherm / p.zdisk, linewidth = 0.9, color = '#F14628')

# ax2 = ax.twinx()
# ax2.plot(d['z'] / p.zdisk, d['compfr'], color = '#A1ABB0', label = '$\\Lambda_C / (\\Lambda_B + \\Lambda_C)$')
# ax2.set_ylim(0,1)
# ax2.set_yscale('linear')
#
# Li1, La1 = ax.get_legend_handles_labels()
# Li2, La2 = ax2.get_legend_handles_labels()
# ax.legend(Li1 + Li2, La1 + La2, fontsize = 9, loc = 'best')

ax.legend(fontsize = 9, loc = 'best')
setx(ax)


#--------------------------------------------------------------------------#

ax = axes[1,2]
ax.set_title(u'Dimensionless variables')

ax.plot(d['z'] / p.zdisk, d['taues'], color = '#AAB5B0', label = '$\\tau_{{\\rm es}}$')
ax.plot(d['z'] / p.zdisk, d['beta'], color = '#245ECD', label = '$\\beta_{{\\rm mag}}$')
ax.set_yscale('log')
ax.set_ylim(1e4,1e-4)

ax2 = ax.twinx()
# ax2.set_ylim(1.2,1.8)
# ax2.axhline(4.0 / 3, linestyle = '--', linewidth = 0.9, color = '#F1C0B5')
# ax2.axhline(5.0 / 3, linestyle = '--', linewidth = 0.9, color = '#B5F1B5')
# ax2.plot(d['z'] / p.zdisk, d['adiab1'], color = '#EC691D', label = '$\\gamma$')
ax2.plot(d['z'] / p.zdisk, d['compy'], label = u'total Y', color = '#D314E6')
# ax2.set_ylim(0,0.25)

lines, labels = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax.legend(lines + lines2, labels + labels2, fontsize = 9, loc = 'best')

ax.axvline(p.zphot  / p.zdisk,  linewidth = 0.9, color = '#4664D0')
ax.axvline(p.ztherm / p.zdisk, linewidth = 0.9, color = '#F14628')

setx(ax)

#--------------------------------------------------------------------------#

Mx = lambda x: (x[1:] + x[:-1]) / 2
Dx = lambda x: (x[1:] - x[:-1])

ax = axes[0,2]
ax.set_title(u'Density')
ax.plot(Mx(d['z'] / p.zdisk), - 2 * Dx(d['pgas']) / (p.omega**2 * Dx(d['z']**2) * cgs_mhydr), color = '#F27463', linewidth = 1.2)
ax.plot(Mx(d['z'] / p.zdisk), - 2 * Dx(d['prad']) / (p.omega**2 * Dx(d['z']**2) * cgs_mhydr), color = '#AAEC96', linewidth = 1.2)
ax.plot(Mx(d['z'] / p.zdisk), - 2 * Dx(d['pmag']) / (p.omega**2 * Dx(d['z']**2) * cgs_mhydr), color = '#77B8F2', linewidth = 1.2)
ax.plot(d['z'] / p.zdisk, d['rho'] / cgs_mhydr, color = '#1D3608', linewidth = 1.8)
ax.set_yscale('log')
ax.set_ylim(1e-7 * max(d['rho']) / cgs_mhydr, 1.5 * max(d['rho']) / cgs_mhydr)

ax.axvline(p.zphot  / p.zdisk,  linewidth = 0.9, color = '#4664D0')
ax.axvline(p.ztherm / p.zdisk, linewidth = 0.9, color = '#F14628')

setx(ax)

#--------------------------------------------------------------------------#

plt.subplots_adjust(0.05, 0.08, 0.95, 0.90, 0.23, 0.19)
plt.savefig('cute.{}'.format(figext))
