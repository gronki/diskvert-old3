from sympy import *
from sympy.printing import *

var('r,rho,T,H')
var('mbh,mdot,alpha')
var('cgs_rschw_sun,cgs_pi,cgs_mcrit_sun,cgs_graw,cgs_msun,cgs_kapes,kappa_0')
var('cgs_stef,cgs_boltz,cgs_mhydr,cgs_a')

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


A = Matrix(fun(r,rho,T,H))
M = A.jacobian([rho,T,H])

ffcode = lambda l,r: fcode(r, assign_to = l, \
    source_format = 'free', standard = 2008)

print(ffcode(MatrixSymbol('A',3,1), A))
print(ffcode(MatrixSymbol('M',3,3), M))
