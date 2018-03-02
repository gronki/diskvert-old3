from diskvert import *
import numpy as np
import matplotlib.pyplot as plt
from par import *

# mbh = 10
# rschw = 2 * cgs_graw * mbh * cgs_msun / cgs_c**2
# omega = np.sqrt(cgs_graw * mbh * cgs_msun / rschw**3)

#------------------------------------------------------------------------------#

T_rad = 1e6
T0 = np.logspace(np.log10(T_rad) + 0.1, np.log10(T_rad) + 2.8, 100)
rho_eqbc = 1e-34 * T_rad**4.5
rho0 = np.logspace(-3.0, 3.5, 130) * rho_eqbc

rho, T= np.meshgrid(rho0, T0)

kappa_0 = 37 * 6.3e22
kappa_es = 0.34
kff = lambda rho,T: kappa_0 * rho * T**(-3.5)
kff_rho = lambda rho, T: kff(rho,T) / rho
kff_T = lambda rho, T: -3.5 * kff(rho,T) / T

KFF = kff(rho,T) * (T + T_rad) * (T**2 + T_rad**2)
KSC = kappa_es * 4 * cgs_boltz * T_rad**4 / (cgs_mel * cgs_c**2)
L = 4 * cgs_stef * rho * (T - T_rad) * (KFF  + KSC)
dL = 2 * cgs_stef * rho / L * ((9 * T_rad**4 - T**4) * kff(rho,T) \
    + 8 * 0.34 * T_rad**4 * cgs_boltz * T / (cgs_mel * cgs_c**2) ) - 1

#------------------------------------------------------------------------------#

plt.figure(figsize = (10,6))
plt.loglog()
c1 = plt.contourf(rho, T, np.log10(L), 12, cmap = 'viridis')
plt.colorbar(c1)
c2 = plt.contourf(rho, T, dL, colors = '#E4E2D9', levels = [-10, 0.0], alpha = 0.18)
c3 = plt.contour(rho, T, dL, colors = '#E4E2D9', levels = [0.0], linewidths = 1.5, alpha = 0.7)
# plt.contourf(rho, T, KSC / (KFF + KSC), levels = [0.0,0.02, 0.10, 0.50, 0.90, 0.98, 1.0], cmap = 'Reds', alpha = 0.2)
c4 = plt.contour(rho, T, KSC / (KFF + KSC), levels = [0.05, 0.50, 0.95], linewidths = [0.4, 1.0, 1.7], colors = '#E4270D', linestyles = '--')
plt.clabel(c4, inline = 1)
c5 = plt.contour(rho, T, np.log10(2 * rho * T * cgs_boltz / cgs_mhydr), levels = np.arange(5,16), colors = '#081A37', linewidths = np.linspace(0.6,1.6,11))
plt.clabel(c5, inline = 1, fmt = '$10^{{%.0f}}$')
# for i in range(4):
#     d,p = col2python('data/A{}.dat'.format(i+1))
#     plt.plot(d['rho'], d['temp'], linewidth = 2)

a = 8.0 / 3.0 * (3.0 / 11.0)**0.125 * kappa_es / kappa_0 * cgs_boltz / (cgs_mel * cgs_c**2)
print a, 5.576e-35
rho_3 = 5.576e-35 * T_rad**4.5

rho1 = np.logspace(-1.2,-0.3) * rho_3
plt.plot(rho1, (8 * kappa_es * cgs_boltz * T_rad**5 )**2 / (3 * kappa_0 * rho1 * cgs_mel * cgs_c**2)**2, color = '#D558EE', linewidth = 2.3, linestyle = '--')
print (8 * kappa_es * cgs_boltz )**2 / (3 * kappa_0 * cgs_mel * cgs_c**2)**2 * 1e60

rho1 = np.logspace(0.5,2.5) * rho_3
plt.plot(rho1, rho1 * 0 + (11.0 / 3)**0.25 * T_rad, color = '#D558EE', linewidth = 2.3, linestyle = '--')
print (11.0 / 3)**0.25

plt.subplots_adjust(right = 0.99)
plt.savefig('cooling.{}'.format(figext))

#------------------------------------------------------------------------------#


# plt.show()
