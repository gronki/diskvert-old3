from sympy import *

var('cgs_c cgs_kapes cgs_k_over_mec2 cgs_stef')
var('rho tgas trad heat')

kabs = Function('kabs')(rho,tgas)
ksct = Function('ksct')(rho,tgas)

yyb = kabs * (trad + tgas) * (trad**2 + tgas**2)
yyc = 4 * ksct * cgs_k_over_mec2 * trad**4
yy = 4 * cgs_stef * rho * (tgas - trad) * (yyb + yyc)

kabsv = IndexedBase('kabsv')
ksctv = IndexedBase('ksctv')

def ff(x,y):
    x = x.subs(Derivative(ksct, rho),  ksctv[2])
    x = x.subs(Derivative(ksct, tgas), ksctv[3])
    x = x.subs(ksct, ksctv[1])
    x = x.subs(Derivative(kabs, rho),  kabsv[2])
    x = x.subs(Derivative(kabs, tgas), kabsv[3])
    x = x.subs(kabs, kabsv[1])
    return fcode(simplify(x), assign_to = y, source_format = 'free',
            standard = 2008, contract = False)

print ff(yy,            'cool')
print ff(yy.diff(rho),  'cool_dr')
print ff(yy.diff(tgas), 'cool_dT')
