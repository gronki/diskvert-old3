from diskvert.cgs import *
import numpy as np

mbh = 1e8
mdot = 0.4
alpha = 0.2

cgs_rschw_sun = 2 * cgs_graw * cgs_msun / cgs_c**2
cgs_mcrit_sun = 4 * cgs_pi * (cgs_graw * cgs_msun / cgs_c) \
        * (cgs_mhydr / cgs_thomson)

r12 = 18 * (alpha * mbh)**(2.0/21) * mdot**(16.0/21)
r23 = 2.5e3 * mdot**(2.0/3)

kappa_0 = 6.13e22

kappa = lambda rho,T: cgs_kapes + kappa_0 * rho * T**(-3.5)

def fun(r,rho,T,H):
    r_cgs = r * mbh * cgs_rschw_sun
    mdot_cgs = mdot * mbh * cgs_mcrit_sun
    omega = ( cgs_graw * cgs_msun * mbh / r_cgs**3 )**0.5
    f = 1 - (3 / r)**0.5
    A1 = 4 * cgs_pi * alpha * rho * H**3 * omega - mdot_cgs * f
    A2 = 0.5625 * alpha * rho**2 * H**4 * omega**3 * kappa(rho,T) - 4 * cgs_stef * T**4
    A3 = 2 * cgs_boltz * rho * T / cgs_mhydr + cgs_a / 3 * T**4 - rho * (omega * H)**2
    return (A1,A2,A3)
