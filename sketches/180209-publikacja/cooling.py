from diskvert import *
import numpy as np
import matplotlib.pyplot as plt
from par import *


#------------------------------------------------------------------------------#

T_rad = 1e6
kappa_0 = 37 * 6.3e22
kappa_es = 0.34
manual_labels = False

T0 = np.logspace(np.log10(T_rad) + 0.01, np.log10(T_rad) + 3.0, 100)
rho_eqbc = 1e-34 * T_rad**4.5
rho0 = np.logspace(-3.0, 3.0, 130) * rho_eqbc

#------------------------------------------------------------------------------#

rho, T = np.meshgrid(rho0, T0)

kff = lambda rho,T: kappa_0 * rho * T**(-3.5)
kff_rho = lambda rho, T: kff(rho,T) / rho
kff_T = lambda rho, T: -3.5 * kff(rho,T) / T

KFF = kff(rho,T) * (T + T_rad) * (T**2 + T_rad**2)
KSC = kappa_es * 4 * cgs_boltz * T_rad**4 / (cgs_mel * cgs_c**2)
L = 4 * cgs_stef * rho * (T - T_rad) * (KFF  + KSC)
dL = 2 * cgs_stef * rho / L * ((9 * T_rad**4 - T**4) * kff(rho,T) \
    + 8 * kappa_es * T_rad**4 * cgs_boltz * T / (cgs_mel * cgs_c**2) ) - 1

#------------------------------------------------------------------------------#

plt.figure(figsize = (figw_sm, figw_sm * 0.77))
plt.loglog()
plt.title(u'The Net Cooling Function ($\\log \\Lambda_{{\\rm rad}}$)')
plt.xlabel(u'$\\rho \\  [{{\\rm g}} / {{\\rm cm}}^3]$')
plt.ylabel(u'$T \\  [{{\\rm K}}]$')

cm = plt.get_cmap('Spectral_r')
c1 = plt.contourf(rho, T, np.log10(L), 13, cmap = cm)
c1b = plt.contour(rho, T, np.log10(L), c1.levels, colors = [cm(0.5)], alpha = 0.5, linewidths = 1.0)

# cb = plt.colorbar(c1, label = '$\\log \\Lambda_{{\\rm rad}}^{{\\rm net}}$')
# cb.add_lines(c1b)


#------------------------------------------------------------------------------#
# compton / brehmstrahlung ratio

fadealpha = lambda color, alphas: ['{}{:02x}'.format(color, int(round(x * 255))) for x in alphas]

# c4 = plt.contourf(rho, T, KSC / (KFF + KSC), levels = [0.50, 1.0], colors = '#FF3D0028')
c4 = plt.contour(rho, T, KSC / (KFF + KSC), levels = [0.05, 0.50, 0.95], linewidths = [0.8, 1.3, 1.8], colors = '#F53B0A', linestyles = '--')
plt.clabel(c4, inline = 1, fmt = '%.2f', manual = manual_labels, colors = '#2F0A00')
# '$\\Lambda_C / \\Lambda$ = %.2f'

num_pressure_contours = 8
c5 = plt.contour(rho, T, np.log10(2 * rho * T * cgs_boltz / cgs_mhydr), num_pressure_contours, colors = '#081A37', linewidths = np.linspace(0.6, 1.0, num_pressure_contours), linestyles = '-.')
# plt.clabel(c5, inline = 1, fmt = '$10^{{%.0f}}$', fontsize = 8, manual = manual_labels)

#------------------------------------------------------------------------------#
# dL / dT

c2 = plt.contourf(rho, T, dL, colors = '#FFFFFF50', levels = [-10, 0.0])
c3 = plt.contour(rho, T, dL, colors = 'white', levels = [0.0], linewidths = 3.0, linestyles = '-')
# plt.clabel(c3, [0.0], fmt = {0.0: '$d \\ln \\Lambda \\ / \\ d \\ln T = 0$'}, inline = 0, manual = manual_labels)

#------------------------------------------------------------------------------#
# rysujemy 5 modeli
nn5 = len(mdots_5)
dp = [ col2python('data/05-' + ix2fn('M', (i,)) + '.dat') for i in range(nn5) ]
rho_ref_0 = 5.576e-35 * T_rad**4.5
cm = plt.get_cmap('Purples')
for i,dp_ in enumerate(dp):
    d,p = dp_
    rho_ref = 5.576e-35 * d['trad']**4.5
    cs = plt.plot(d['rho'] * rho_ref_0 / rho_ref, d['temp'] * T_rad / d['trad'], linewidth = (0.7 + 1.1 * i / (nn5 - 1.0)), color = '#0E0FBB', linestyle = '-')

#------------------------------------------------------------------------------#

a1 = (11.0 / 3)**0.25
a2 = (8 * kappa_es * cgs_boltz)**2 / (3 * kappa_0 * cgs_mel * cgs_c**2)**2 * (1e6)**10
a3 = 8.0 / 3.0 * (3.0 / 11.0)**0.125 * kappa_es / kappa_0 * cgs_boltz / (cgs_mel * cgs_c**2) * (1e6)**4.5

a = 8.0 / 3.0 * (3.0 / 11.0)**0.125 * kappa_es / kappa_0 * cgs_boltz / (cgs_mel * cgs_c**2)
print a, 5.576e-35, a3
rho_3 = 5.576e-35 * T_rad**4.5

rho1 = np.logspace(-3,0) * rho_3
# plt.plot(rho1, (8 * kappa_es * cgs_boltz * T_rad**5 )**2 / (3 * kappa_0 * rho1 * cgs_mel * cgs_c**2)**2, color = '#383838', linewidth = 1.6, linestyle = ':')
print (8 * kappa_es * cgs_boltz )**2 / (3 * kappa_0 * cgs_mel * cgs_c**2)**2 * 1e60, a2

rho1 = np.logspace(0,4) * rho_3
# plt.plot(rho1, rho1 * 0 + (11.0 / 3)**0.25 * T_rad, color = '#383838', linewidth = 1.6, linestyle = ':')
print (11.0 / 3)**0.25, a1

#------------------------------------------------------------------------------#

whitebox = dict(boxstyle = 'square', fc = '#FFFFFFA0', ec = 'none')

plt.annotate(u'$d \\ln \\Lambda \\ / \\ d \\ln T > 0$\n(stable)',
        (0.03, 0.24), xycoords = 'axes fraction', color = 'black', fontsize = 11, bbox = whitebox)
plt.annotate(u'$d \\ln \\Lambda \\ / \\ d \\ln T < 0$\n(unstable)',
        (0.44, 0.44), xycoords = 'axes fraction', color = 'black', fontsize = 11, bbox = whitebox)

plt.annotate(u'$\\Lambda_C > \\Lambda_B$', (0.51, 0.70), xycoords = 'axes fraction', color = '#C72C03', fontsize = 11, bbox = whitebox)
plt.annotate(u'$\\Lambda_C < \\Lambda_B$', (0.71, 0.65), xycoords = 'axes fraction', color = '#C72C03', fontsize = 11, bbox = whitebox)

#------------------------------------------------------------------------------#

plt.xlim(np.min(rho), np.max(rho))
plt.ylim(np.min(T), np.max(T))

plt.subplots_adjust(left = 0.11, bottom = 0.13, right = 0.96, top = 0.92)
plt.savefig('cooling.{}'.format(figext))

#------------------------------------------------------------------------------#


# plt.show()
