# coding: utf-8
from sympy import *
var('sigma kappa_0 T rho T_rad kappa_es k m_el c tmp kappa_ff', real = True, positive = True)

fix = lambda x: x.subs(kappa_0, kappa_ff / (rho * T**Rational(-7,2)))

lam1 = 4 * sigma * kappa_0 * rho**2 * T**Rational(-7,2) * (T**4 -  T_rad**4)
lam1_dT = lam1.subs(rho, tmp / T).diff(T).subs(tmp, rho * T).simplify()
print 'Lff = \n', pretty(fix(lam1))
print 'dLff / dT = \n', pretty(fix(lam1_dT))
Tcrit_1 = simplify(solve(lam1_dT, T))
print 'dLff / dT = 0  when  T = \n', pretty(Tcrit_1)

lam2 = 16 * sigma * T_rad**4 * kappa_es * rho * k * (T - T_rad) / (m_el * c**2)
lam2_dT = lam2.subs(rho, tmp / T).diff(T).subs(tmp, rho * T).simplify()
print 'Les = \n', pretty(simplify(fix(lam2)))
print 'dLes / dT = \n',  pretty(simplify(fix(lam2_dT)))

lam = lam1 + lam2
lam_dT = lam1_dT + lam2_dT
print 'Ltot = \n', pretty(simplify(fix(lam)))
print 'dLtot / dT = \n',  pretty(simplify(fix(lam_dT)))
print lam_dT.subs(T_rad,0)
Tcrit_2 = solve(lam1_dT.subs(T_rad,0) + lam2_dT, T)
print 'dLtot / dT = 0  when  T = \n',  pretty(Tcrit_2)

rhocrit_12 = simplify(solve(Tcrit_1[0] - Tcrit_2[0], rho))
print 'Tcrit_1 = Tcrit_2 for rho = \n', (rhocrit_12)
