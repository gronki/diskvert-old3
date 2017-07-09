from sympy import symbols, Function, pretty_print, IndexedBase, Symbol
from sympy.printing.fcode import FCodePrinter

fcprinter = FCodePrinter(dict(
    source_format = 'free',
    standard = 2008,
    user_functions = {
        'absiso': 'absiso',
        'sctiso': 'sctiso',
    }
))
def fcode(rhs, assign_to = None):
    return fcprinter.doprint(rhs, assign_to)

def f2vec(expr, functions, fval = IndexedBase('F')):
    from sympy import Wild, Derivative
    for ifun,f in zip(range(len(functions)),functions):
        nargs = len(f.args)
        ff = f.func
        w = [ Wild('w{}'.format(i)) for i in range(nargs) ]
        for i in range(nargs):
            replace_what = Derivative(ff(*w),w[i])
            expr = expr.replace( replace_what, fval[i+1,ifun] )
        expr = expr.replace( ff(*w), fval[0,ifun] )
    return expr.doit()

T, Trad, Pgas = symbols('T Trad Pgas', real = True, positive = True)
heat = symbols('heat', real = True, positive = True)
sigma, c, m_el, k, mH, miu = symbols('cgs_stef cgs_c cgs_mel cgs_boltz cgs_mhydr miu', real = True, positive = True)

rho = Symbol('rhotemp') / T

kappa_abs = Function('absiso')(T)
kappa_sct = Function('sctiso')(T)

cool = 4 * sigma * rho * ( T - Trad ) * (
            kappa_abs * (T + Trad) * (T**2 + Trad**2) \
            + kappa_sct * 4 * k * Trad**4 / ( m_el * c**2 )
        )

bil = cool - heat
bil = bil.subs(k / ( m_el * c**2 ), Symbol('cgs_k_over_mec2'))

from sympy import MatrixSymbol
fs = [kappa_abs,kappa_sct]
fval = MatrixSymbol('F', 2, len(fs))

print fcode(f2vec(bil.simplify(), fs, fval), 'bil')
print fcode(f2vec(bil.diff(T).simplify(), fs, fval), 'bil_T')
