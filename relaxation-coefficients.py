#!/usr/bin/env python
# coding: utf-8

from sympy import Symbol, Function, Lambda, Derivative, IndexedBase, Eq, symbols, Integer, Rational, Matrix, MatrixSymbol, Wild
from numpy import ndarray, zeros, linspace, logspace, meshgrid
from sys import stdin, stdout, stderr

#------------------------------------------------------------------------------#

alpha   = Symbol('alpha', real = True, positive = True)
eta     = Symbol('eta', real = True, positive = True)
omega   = Symbol('omega', real = True, positive = True)
F_acc   = Symbol('F_acc', real = True, positive = True)
mu      = Symbol('mu', real = True, positive = True)
m_el     = Symbol('m_el', real = True, positive = True)
m_H      = Symbol('m_H', real = True, positive = True)
c        = Symbol('c', real = True, positive = True)
sigma    = Symbol('sigma', real = True, positive = True)
k        = Symbol('k_B', real = True, positive = True)

global_variables = {
    omega:  'omega',
    F_acc:  'facc',
    m_el:   'cgs_mel',
    m_H:    'cgs_mhydr',
    c:      'cgs_c',
    sigma:  'cgs_stef',
    k:      'cgs_boltz',
    mu:     'miu',
    alpha:  'alpha',
    eta:    'zeta',
}

#------------------------------------------------------------------------------#

from sympy.printing.fcode import FCodePrinter
from sympy.printing.precedence import precedence, PRECEDENCE

class F90CodePrinter(FCodePrinter):
    def __init__(self, settings = {}):
        settings['source_format'] = 'free'
        settings['standard'] = 2008
        super(F90CodePrinter, self).__init__(settings)
        self._relationals['!='] = '.ne.'
        self._relationals['=='] = '.eq.'
        self._relationals['>='] = '.ge.'
        self._relationals['>'] = '.gt.'
        self._relationals['<='] = '.le.'
        self._relationals['<'] = '.lt.'
        self._relationals['&&'] = '.and.'
        self._relationals['||'] = '.or.'
        self._relationals['!'] = '.not.'
        self._lead_cont = ' '*3 + '&' + ' '*3

    def _print_MatrixElement(self, expr):
        if expr.parent.shape[1] <= 1:
            return "{0}({1})".format(expr.parent, expr.i + 1)
        return "{0}({1},{2})".format(expr.parent, expr.i + 1, expr.j + 1)

fprinter = F90CodePrinter({'source_format':'free', 'standard':2008})

#------------------------------------------------------------------------------#

fcoeffsrc = open('src/coefficients.fi','w')

fsub_coeff = """pure subroutine {name} (z,Y,D,F,A,MY,MD)
use iso_fortran_env, only: real64
implicit none
real(real64), intent(in) :: z
! first column: values, second column: 1ord derivatives
real(real64), dimension(:), intent(in) :: Y,D
! each row is one function: its value and derivatives to rho and T
real(real64), dimension(:,:), intent(in) :: F
! right-hand-side of the equation
real(real64), dimension(:), intent(out) :: A
! jacobians with respect to Y and dY/dz
real(real64), dimension(:,:), intent(out) :: MY, MD
{A}\n{MY}\n{MD}\nend subroutine\n\n"""

fsub_bound = """pure subroutine {name}_{lr} (z,Y,F,B,MB)
use iso_fortran_env, only: real64
implicit none
real(real64), intent(in) :: z
real(real64), dimension(:), intent(in) :: Y
real(real64), dimension(:,:), intent(in) :: F
real(real64), dimension(:), intent(out) :: B
real(real64), dimension(:,:), intent(out) :: MB
{B}\n{MB}\nend subroutine\n\n"""

#------------------------------------------------------------------------------#

fswptrs = open('src/mrxptrs.fi','w')
fswname = open('src/mrxname.fi','w')
fswdims = open('src/mrxdims.fi','w')
fswhash = open('src/mrxhash.fi','w')
fswynum = open('src/mrxynum.fi','w')
fswall = [fswptrs, fswdims, fswhash, fswname, fswynum]

for f in fswall: f.write("select case (nr)\n")

#------------------------------------------------------------------------------#

for enableCorona, enableMagnetic, enableConduction in [
    (False, False, False),
    (True,  False, False),
    (False, True,  False),
    (True,  True,  False),
    (False, False, True),
    (True,  False, True),
    (False, True,  True),
    (True,  True,  True),
]:

    # współrzędna, po której liczymy
    z       = Symbol('z')

    # niewiadome
    rho     = Function('rho')(z)
    T_gas   = Function('T_gas' if enableCorona else 'T')(z)
    F_rad   = Function('F_rad')(z)

    T_rad   = Function('T_rad')(z) if enableCorona else T_gas
    P_mag   = Function('P_mag')(z) if enableMagnetic else 0
    F_cond  = Function('F_cond')(z) if enableConduction else 0

    #--------------------------------------------------------------------------#

    # niewiadome w postaci listy
    yvar = [ rho, T_gas, F_rad ]
    if enableCorona: yvar.append(T_rad)
    if enableMagnetic: yvar.append(P_mag)
    if enableConduction: yvar.append(F_cond)

    yvar_all = [rho, T_gas, T_rad, F_rad, P_mag, F_cond]
    yval_hash = [ yvar.index(y) if y in yvar else -100 for y in yvar_all ]

    #--------------------------------------------------------------------------#

    # nieprzezroczystosci i rpzewodnizctwa
    kappa_abs    = Function('fkabs')(rho,T_gas)
    kappa_sct    = Function('fksct')(rho,T_gas)
    kappa_cond   = Function('fkcnd')(rho,T_gas)
    global_functions = set([ kappa_abs, kappa_sct, kappa_cond ])

    fval = MatrixSymbol('F', 3, len(global_functions))

    def function_derivatives(eq):
        functions_extra = []
        for ifun,f in zip(range(len(global_functions)),global_functions):
            nargs = len(f.args)
            ff = f.func
            w = [ Wild('w{}'.format(i)) for i in range(nargs) ]
            for i in range(nargs):
                replace_what = Derivative(ff(*w),w[i])
                eq = eq.replace( replace_what, fval[i+1,ifun] )
            eq = eq.replace( ff(*w), fval[0,ifun] )
        return eq.doit()

    #--------------------------------------------------------------------------#

    # calkowita nieprzezroczystosc
    kappa = kappa_abs + kappa_sct
    # predkosc wznoszenia pola
    vrise = Lambda(z, eta * omega * z)
    # cisnienie gazu
    P_gas = k * rho * T_gas / ( mu * m_H )
    # cisnienie promieniowania
    P_rad = 4 * sigma / (3 * c) * T_rad**4
    # cisnienie calkowite
    P_tot = P_gas + P_rad + P_mag
    # strumien Poyntinga
    F_mag = 2 * P_mag * vrise(z)

    #--------------------------------------------------------------------------#

    equations = []
    boundL = []
    boundR = []

    #--------------------------------------------------------------------------#
    # Hydrostatic equilibrium
    #
    equations.append(
        (P_gas.diff(z) - kappa * rho / c * F_rad + Derivative(P_mag,z))
        + omega**2 * rho * z
    )

    #--------------------------------------------------------------------------#
    # Dyfuzja promieniowania
    #
    zeta_rad = 16 * sigma * T_rad**3 / (3 * kappa * rho)
    equations.append(
        F_rad + zeta_rad * Derivative(T_rad,z)
    )


    #--------------------------------------------------------------------------#
    # Bilans grzania i chlodzenia
    #

    GasHeating = alpha * omega * P_tot
    if enableMagnetic: GasHeating = - vrise(z) * Derivative(P_mag,z)
    equations.append(
        Derivative(F_rad,z) + Derivative(F_cond,z) - GasHeating
    )
    boundL.append( F_rad )
    boundR.append( F_rad + F_mag + F_cond - F_acc )
    boundR.append( F_rad - 2 * sigma * T_rad**4 )

    #--------------------------------------------------------------------------#
    # Magnetic field evolution
    #
    if enableMagnetic:
        equations.append(
            2 * P_mag * vrise(z).diff(z) \
            + vrise(z) * Derivative(P_mag,z) \
            - alpha * omega * P_tot
        )
        boundL.append( 2 * P_mag * vrise(z).diff(z) - alpha * omega * P_tot )

    #--------------------------------------------------------------------------#
    # Funkcja zrodlowa
    #
    if enableCorona:
        SourceFunction = 4 * sigma * rho * ( T_gas - T_rad ) * (
            kappa_abs * (T_gas + T_rad) * (T_gas**2 + T_rad**2) \
            + kappa_sct * 4 * k * T_rad**4 / ( m_el * c**2 )
        )
        equations.append(SourceFunction - Derivative(F_rad,z))
        boundL.append( T_gas - T_rad )

    #--------------------------------------------------------------------------#
    # Dyfuzja ciepla przez przewodnictwo
    #
    if enableConduction:
        equations.append(
            F_cond + kappa_cond * Derivative(T_gas,z)
        )
        boundL.append( F_cond )

    #--------------------------------------------------------------------------#

    assert len(equations) == len(yvar)
    assert len(boundL) + len(boundR) == len(yvar)

    neq = len(equations)
    ny = len(yvar)
    nbl = len(boundL)
    nbr = len(boundR)

    Y = MatrixSymbol('Y', ny, 1)
    D = MatrixSymbol('D', ny, 1)

    #--------------------------------------------------------------------------#

    def discretize(eq):
        if eq == 0: return eq
        for iy,y in zip(range(len(yvar)),yvar):
            eq = eq.subs(Derivative(y,z),D[iy]).subs(y,Y[iy])
        for k,v in global_variables.items():
            eq = eq.subs(k,Symbol(v, real = True, positive = k.is_positive))
        return function_derivatives(eq)

    #--------------------------------------------------------------------------#

    basename = '{magn}{comp}{cond}'.format(
        comp = 'cmpt' if enableCorona else 'dyfu',
        magn = 'magn' if enableMagnetic else 'alph',
        cond = 'cond' if enableConduction else '',
    )
    routine_name = "mrx_coeff" + basename
    model_nr = 1 \
        + (1 if enableCorona else 0) \
        + (2 if enableMagnetic else 0) \
        + (4 if enableConduction else 0)

    #--------------------------------------------------------------------------#

    A = Matrix(zeros((neq,)))
    MY = Matrix(zeros((neq,ny)))
    MD = Matrix(zeros((neq,ny)))

    for ieq,eq in zip(range(len(equations)),equations):
        A[ieq] = discretize(-eq)
        for iy,y in zip(range(len(yvar)),yvar):
            MY[ieq,iy] = discretize(eq.diff(y))
            MD[ieq,iy] = discretize(eq.diff(Derivative(y,z)))

    #--------------------------------------------------------------------------#

    BL = Matrix(zeros((nbl,)))
    MBL = Matrix(zeros((nbl,ny)))

    for ib,b in zip(range(len(boundL)),boundL):
        BL[ib] = discretize(-b)
        for iy,y in zip(range(len(yvar)),yvar):
            MBL[ib,iy] = discretize(b.diff(y))

    #--------------------------------------------------------------------------#

    BR = Matrix(zeros((nbr,)))
    MBR = Matrix(zeros((nbr,ny)))

    for ib,b in zip(range(len(boundR)),boundR):
        BR[ib] = discretize(-b)
        for iy,y in zip(range(len(yvar)),yvar):
            MBR[ib,iy] = discretize(b.diff(y))

    #--------------------------------------------------------------------------#

    fcoeffsrc.write("\n!{s}!\n! {A}\n! {B}\n! {C}\n!{s}!\n\n".format(
        s = "-"*78,
        A = ("heating-cooling balance" if enableCorona else "thermal diffusion"),
        B = ("magnetic heating" if enableMagnetic else "alpha prescription"),
        C = ("radiation + thermal conduction" if enableConduction \
            else "radiative energy transport"),
    ))

    fcoeffsrc.write(fsub_coeff.format(
        name = routine_name,
        ny = ':', na = ':', nf = len(global_functions),
        A = fprinter.doprint(A, 'A'),
        MY = fprinter.doprint(MY, 'MY'),
        MD = fprinter.doprint(MD, 'MD'),
    ))
    fcoeffsrc.write(fsub_bound.format(
        name = routine_name, lr = 'bl',
        ny = ':', nb = ':', nf = len(global_functions),
        B = fprinter.doprint(BL, 'B'),
        MB = fprinter.doprint(MBL, 'MB'),
    ))
    fcoeffsrc.write(fsub_bound.format(
        name = routine_name, lr = 'br',
        ny = ":", nb = ":", nf = len(global_functions),
        B = fprinter.doprint(BR, 'B'),
        MB = fprinter.doprint(MBR, 'MB'),
    ))

    fcoeffsrc.write("\n")

    #--------------------------------------------------------------------------#

    for f in fswall: f.write('case({})\n'.format(model_nr))

    fswname.write ('    name = \"{}\"\n'.format(basename))
    fswptrs.write ('    f => {}\n'.format(routine_name))
    fswptrs.write ('    fbL => {}\n'.format(routine_name+'_bL'))
    fswptrs.write ('    fbR => {}\n'.format(routine_name+'_bR'))
    fswynum.write ('    ny = {}\n'.format(len(yvar)))
    fswdims.write ('    n = {}\n'.format(len(yvar)))
    fswdims.write ('    nbL = {}\n'.format(len(boundL)))
    fswdims.write ('    nbR = {}\n'.format(len(boundR)))
    fswhash.write ('    ihash = [{}]\n'.format(",".join([ str(i+1) \
            for i in yval_hash ])))

#------------------------------------------------------------------------------#

for f in fswall:
    f.write("case default\n    error stop\nend select\n")
    f.close()

fcoeffsrc.close()
