from sympy import *

var('cgs_c cgs_kapes cgs_k_over_mec2 cgs_stef')
var('rho tgas trad heat')

kabp = Function('kabp')(rho,tgas)
ksct = Function('ksct')(rho,tgas)

yyb = kabp * (tgas**4 - trad**4)
# tcompt = sqrt((tgas)**2 + (4 * cgs_k_over_mec2 * tgas**2)**2)
yyc = ksct * trad**4 * cgs_k_over_mec2 * 4 * (tgas - trad)
yy = 4 * cgs_stef * rho * (yyb + yyc)

kabpv = IndexedBase('kabpv')
ksctv = IndexedBase('ksctv')

def ff(x,y):
    x = x.subs(Derivative(ksct, rho),  ksctv[2])
    x = x.subs(Derivative(ksct, tgas), ksctv[3])
    x = x.subs(ksct, ksctv[1])
    x = x.subs(Derivative(kabp, rho),  kabpv[2])
    x = x.subs(Derivative(kabp, tgas), kabpv[3])
    x = x.subs(kabp, kabpv[1])
    return fcode(simplify(x), assign_to = y, source_format = 'free',
            standard = 2008, contract = False)

print ff(yy,            'cool')
print ff(yy.diff(rho),  'cool_dr')
print ff(yy.diff(tgas), 'cool_dT')

def gg(x):
    x = x.subs(Derivative(ksct, rho),  0)
    x = x.subs(Derivative(ksct, tgas), 0)
    x = x.subs(Derivative(kabp, rho),  kabp / rho)
    x = x.subs(Derivative(kabp, tgas), - Rational(7,2) * kabp / tgas)
    return simplify(x)

print gg(yy.diff(tgas) - rho / tgas * yy.diff(rho))
