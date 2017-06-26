from diskvert.cgs import *
import numpy as np

mbh = 1e8
mdot = 0.1
alpha = 0.3

kappa_ff = 6.13e22
kappa = lambda rho,T: cgs_kapes + kappa_ff * rho * T**(-3.5)
P = lambda rho,T: 2 * cgs_boltz * rho * T / cgs_mhydr + cgs_a / 3 * T**4

def fun(mbh,mdot,r,alpha,rho,T,H):
    r_cgs = r * mbh * sol_rschw
    mdot_cgs = mdot * mbh * sol_mdot_edd
    omega = ( cgs_graw * cgs_msun * mbh / r_cgs**3 )**0.5
    f = 1 - (3 / r)**0.5
    A1 = 4 * cgs_pi * alpha * rho * H**3 * omega - mdot_cgs * f
    A2 = 0.5625 * alpha * (rho * H**2)**2 * omega**3 * kappa(rho,T) \
            - 4 * cgs_stef * T**4
    A3 = P(rho,T) - rho * (omega * H)**2
    return (A1,A2,A3)

def calc_r123_old(mbh,mdot,alpha):
    r12 = 18 * (alpha * mbh)**(2.0/21) * mdot**(16.0/21)
    r23 = 2.5e3 * mdot**(2.0/3)
    return r12,r23

def calc_rTH_old(mbh,mdot,r,alpha):
    f = 1 - np.sqrt(3 / r)

    rho1 = 9e-4 * (alpha * mbh)**(-1) * mdot**(-2) * r**(3./2) * f**(-2)
    T1  = 4.9e7 * (alpha * mbh)**(-0.25) * r**(-3.0/8)
    H1 = 5.5e4 * mbh * mdot * f + r*0

    rho2 = 8 * (alpha*mbh)**(-0.7) * mdot**(0.4) * r**(-33./20) * f**(0.4)
    T2  = 2.2e8 * (alpha*mbh)**(-0.2) * mdot**(0.4) * r**(-0.9) * f**(0.4)
    H2 = 2.7e3 * alpha**(-0.1) * mbh**(0.9) * mdot**(0.2) * r**1.05 * f**0.2

    rho3 = 47 * (alpha*mbh)**(-0.7) * mdot**0.55 * r**(-1.875) * f**0.55
    T3  = 6.9e7 * (alpha*mbh)**(-0.2) * mdot**(0.3) * r**(-0.75) * f**(0.3)
    H3 = 1.5e3 * alpha**(-0.1) * mbh**0.9 * mdot**0.15 * r**1.125 * f**0.15

    return [ (rho1, rho2, rho3), (T1, T2, T3), (H1, H2, H3) ]

def calc_r123_new(mbh,mdot,alpha):
    r12 = 164.17*alpha**0.095238*cgs_kapes**0.85714*mbh**0.095238*mdot**0.7619*1**0.7619
    r23 = 5.0556e+19*cgs_kapes**1.3333*kappa_ff**(-0.66667)*mdot**0.66667*1**0.66667
    f = lambda r: 1 - sqrt(3/r) if r > 3 else 0
    for i in range(4):
        r12 = 164.17*alpha**0.095238*cgs_kapes**0.85714*mbh**0.095238*mdot**0.7619*f(r12)**0.7619
        r23 = 5.0556e+19*cgs_kapes**1.3333*kappa_ff**(-0.66667)*mdot**0.66667*f(r23)**0.66667
    return r12, r23

def calc_rTH_new(mbh,mdot,r,alpha):
    f = lambda r: np.where(r > 3, 1 - np.sqrt(3/r), 0)

    rho1 = 2.3049e-6*mbh**(-1.0)*r**1.5/(alpha*cgs_kapes**3*mdot**2*f(r)**2)
    T1 = 4.6189e+7*alpha**(-0.25)*cgs_kapes**(-0.25)*mbh**(-0.25)*r**(-0.375)
    H1 = 9.8298e+5*cgs_kapes*mbh**1.0*mdot*f(r)
    rho2 = 21.921*alpha**(-0.7)*cgs_kapes**(-0.3)*mbh**(-0.7)*mdot**0.4*r**(-1.65)*f(r)**0.4
    T2 = 6.7232e+8*alpha**(-0.2)*cgs_kapes**0.2*mbh**(-0.2)*mdot**0.4*r**(-0.9)*f(r)**0.4
    H2 = 4639.5*alpha**(-0.1)*cgs_kapes**0.1*mbh**0.9*mdot**0.2*r**1.05*f(r)**0.2
    rho3 = 5.9458e+5*alpha**(-0.7)*kappa_ff**(-0.15)*mbh**(-0.7)*mdot**0.55*r**(-1.875)*f(r)**0.55
    T3 = 7.4475e+5*alpha**(-0.2)*kappa_ff**0.1*mbh**(-0.2)*mdot**0.3*r**(-0.75)*f(r)**0.3
    H3 = 154.42*alpha**(-0.1)*kappa_ff**0.05*mbh**0.9*mdot**0.15*r**1.125*f(r)**0.15

    return [ (rho1, rho2, rho3), (T1, T2, T3), (H1, H2, H3) ]
