from sympy import *
from parameters import fun
from subprocess import call
from diskvert.factory import *
from os import path

def FF():
    r,rho,T,H = symbols('r,rho,T,H')
    mbh,mdot,alpha = symbols('mbh,mdot,alpha')

    A = Matrix(fun(mbh,mdot,r,alpha,rho,T,H))
    M = A.jacobian([rho,T,H])

    return A,M

def regenerate_coefficient_routine():

    A,M = FF()

    ffcode = lambda l,r: fcode(r, assign_to = l, \
        source_format = 'free', standard = 2008)

    srcfn = 'coeff.f90'
    libfn = 'coeff.so'

    with open(srcfn,'w') as f:
        f.write("subroutine coeffs(mbh,mdot,r,alpha,rho,T,H,A,M) bind(C)\n" \
            "implicit none\n" \
            "double precision, intent(in), value :: mbh,mdot,r,alpha,rho,T,H\n" \
            "double precision, intent(out) :: A(3,1), M(3,3)\n"
            "{A}\n{M}\nend subroutine\n".format(
                A = ffcode(MatrixSymbol('A',3,1), A.evalf(7)),
                M = ffcode(MatrixSymbol('M',3,3), M.evalf(7)),
        ))

    call("f95 -g -pipe -O3 -mieee-fp -shared -fPIC {s} -o {l}".format(
        s = srcfn, l = libfn,
    ).split())

    return CDLL(libfn)

__lib = regenerate_coefficient_routine()

compcoeffs = CFACTORY(__lib, 'coeffs', [
    (c_double, 'in', 'mbh'),
    (c_double, 'in', 'mdot'),
    (c_double, 'in', 'r'),
    (c_double, 'in', 'alpha'),
    (c_double, 'in', 'rho'),
    (c_double, 'in', 'T'),
    (c_double, 'in', 'H'),
    (ndpointer(c_double,1), 'inout', 'A'),
    (ndpointer(c_double,2), 'inout', 'M'),
])
