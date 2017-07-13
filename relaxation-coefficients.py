#!/usr/bin/env python
# coding: utf-8

from sympy import Symbol, Function, Lambda, Derivative, IndexedBase, Eq, symbols, Integer, Rational, Matrix, MatrixSymbol, Wild, simplify
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
use iso_fortran_env, only: dp => real64
implicit none
real(dp), intent(in) :: z
! first column: values, second column: 1ord derivatives
real(dp), dimension(:), intent(in) :: Y,D
! each row is one function: its value and derivatives to rho and T
real(dp), dimension(:,:), intent(in) :: F
! right-hand-side of the equation
real(dp), dimension(:), intent(out) :: A
! jacobians with respect to Y and dY/dz
real(dp), dimension(:,:), intent(out) :: MY, MD
{A}\n{MY}\n{MD}\nend subroutine\n\n"""

fsub_bound = """pure subroutine {name}_{lr} (z,Y,F,B,MB)
use iso_fortran_env, only: dp => real64
implicit none
real(dp), intent(in) :: z
real(dp), dimension(:), intent(in) :: Y
real(dp), dimension(:,:), intent(in) :: F
real(dp), dimension(:), intent(out) :: B
real(dp), dimension(:,:), intent(out) :: MB
{B}\n{MB}\nend subroutine\n\n"""

fsub_constr = """pure subroutine {name}_c (z,Y,F,C,MC)
use iso_fortran_env, only: dp => real64
implicit none
real(dp), intent(in) :: z
real(dp), dimension(:), intent(in) :: Y
real(dp), dimension(:,:), intent(in) :: F
real(dp), dimension(:), intent(out) :: C
real(dp), dimension(:,:), intent(out) :: MC
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

choices = [ (bil,bilf,mag,cnd)          \
            for cnd in [False,True]     \
            for mag in [False,True]     \
            for bil in [False,True]     \
            for bilf in [True,False]    ]

for balance, bilfull, magnetic, conduction in choices:

    if bilfull and not balance: continue

    # współrzędna, po której liczymy
    z       = Symbol('z')

    # niewiadome
    rho     = Function('rho')(z)
    T_gas   = Function('T_gas' if balance else 'T')(z)
    F_rad   = Function('F_rad')(z)

    T_rad   = Function('T_rad')(z) if balance else T_gas
    P_mag   = Function('P_mag')(z) if magnetic else 0
    F_cond  = Function('F_cond')(z) if conduction else 0

    #--------------------------------------------------------------------------#

    # niewiadome w postaci listy
    yvar = [ rho, T_gas, F_rad ]
    if balance: yvar.append(T_rad)
    if magnetic: yvar.append(P_mag)
    if conduction: yvar.append(F_cond)

    yvar_all = [rho, T_gas, T_rad, F_rad, P_mag, F_cond]
    yval_hash = [ yvar.index(y) if y in yvar else -1 for y in yvar_all ]

    #--------------------------------------------------------------------------#

    # nieprzezroczystosci i rpzewodnizctwa
    kabs = Function('fkabs')(rho,T_gas)
    ksct = Function('fksct')(rho,T_gas)
    kcnd = Function('fkcnd')(rho,T_gas)

    fvec = [ kabs, ksct, kcnd ]

    fval = MatrixSymbol('F', 3, len(fvec))

    def function_derivatives(eq):
        functions_extra = []
        for ifun,f in zip(range(len(fvec)),fvec):
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
    kappa = kabs + ksct
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
    constraints = []
    boundL = []
    boundR = []

    #--------------------------------------------------------------------------#
    # Hydrostatic equilibrium
    #
    equations.append(P_gas.diff(z)      \
            - kappa * rho / c * F_rad   \
            + Derivative(P_mag,z)       \
            + omega**2 * rho * z)
    # equations.append(P_tot.diff(z) + omega**2 * rho * z)

    #--------------------------------------------------------------------------#
    # Dyfuzja promieniowania
    #
    equations.append(
        3 * kappa * rho * F_rad + 16 * sigma * T_rad**3 * Derivative(T_rad,z)
    )

    #--------------------------------------------------------------------------#
    # Bilans grzania i chlodzenia
    #
    if magnetic:
        heat = 2 * P_mag * vrise(z).diff(z) - alpha * omega * P_tot
    else:
        heat = alpha * omega * P_tot

    equations.append(
        heat - Derivative(F_rad,z) - Derivative(F_cond,z)
    )

    boundL.append( F_rad )
    boundR.append( F_rad + F_mag + F_cond - F_acc )
    boundR.append( F_rad - 2 * sigma * T_rad**4 )

    #--------------------------------------------------------------------------#
    # Magnetic field evolution
    #
    if magnetic:
        equations.append(
            2 * P_mag * vrise(z).diff(z) \
            + vrise(z) * Derivative(P_mag,z) \
            - alpha * omega * P_tot
        )
        boundL.append( 2 * P_mag * vrise(z).diff(z) - alpha * omega * P_tot )

    #--------------------------------------------------------------------------#
    # Funkcja zrodlowa
    #
    if balance:
        sdyf = kabs * (T_gas + T_rad) * (T_gas**2 + T_rad**2)
        ssct = ksct * 4 * k * T_rad**4 / ( m_el * c**2 )
        radcool = 4 * sigma * rho * (T_gas - T_rad) * \
            ((sdyf if bilfull else 0) + ssct)
        if conduction:
            equations.append(radcool - Derivative(F_rad,z))
            boundL.append(T_gas - T_rad)
        else:
            constraints.append(heat - radcool)

    #--------------------------------------------------------------------------#
    # Dyfuzja ciepla przez przewodnictwo
    #
    if conduction:
        equations.append(
            F_cond + kcnd * Derivative(T_gas,z)
        )
        boundL.append(F_cond)

    #--------------------------------------------------------------------------#

    na = len(equations)
    nc = len(constraints)
    ny = len(yvar)
    nbl = len(boundL)
    nbr = len(boundR)

    assert na + nc == ny
    assert na == nbl + nbr

    Y = MatrixSymbol('Y', ny, 1)
    D = MatrixSymbol('D', ny, 1)

    #--------------------------------------------------------------------------#

    def discretize(eq):
        if eq == 0: return eq
        for iy,y in zip(range(ny),yvar):
            eq = eq.subs(Derivative(y,z),D[iy]).subs(y,Y[iy])
        for k,v in global_variables.items():
            eq = eq.subs(k,Symbol(v, real = True, positive = k.is_positive))
        return simplify(function_derivatives(eq))

    #--------------------------------------------------------------------------#

    basename = '{magn}{comp}{cond}'.format(
        magn = 'm' if magnetic else 'a',
        comp = ('c' if bilfull else 'w') if balance else  'd',
        cond = 't' if conduction else '',
    )
    model_nr = 1 \
        + (1 if balance else 0) \
        + (2 if bilfull else 0) \
        + (4 if magnetic else 0) \
        + (8 if conduction else 0)

    routine_name = "mrx_coeff_{}".format(basename)
    #--------------------------------------------------------------------------#

    fcoeffsrc.write("! {}\n".format("heating-cooling balance ({})".format("full equation" if bilfull else "compton term only") if balance else "thermal diffusion"))
    fcoeffsrc.write("! {}\n".format("magnetic heating" if magnetic else "alpha prescription"))
    fcoeffsrc.write("! {}\n".format("radiation + thermal conduction" if conduction \
        else "radiative energy transport"))
    #--------------------------------------------------------------------------#

    A = Matrix(zeros((na,)))
    MY = Matrix(zeros((na,ny)))
    MD = Matrix(zeros((na,ny)))

    for ieq,eq in zip(range(na),equations):
        A[ieq] = discretize(-eq)
        for iy,y in zip(range(ny),yvar):
            MY[ieq,iy] = discretize(eq.diff(y))
            MD[ieq,iy] = discretize(eq.diff(Derivative(y,z)))

    fcoeffsrc.write(fsub_coeff.format(
        name = routine_name,
        A = fprinter.doprint(A, 'A'),
        MY = fprinter.doprint(MY, 'MY'),
        MD = fprinter.doprint(MD, 'MD'),
    ))

    #--------------------------------------------------------------------------#

    if nc > 0:
        C = Matrix(zeros((nc,)))
        MC = Matrix(zeros((nc,ny)))

        for ict,ct in zip(range(nc),constraints):
            C[ict] = discretize(-ct)
            for iy,y in zip(range(ny),yvar):
                MC[ict,iy] = discretize(ct.diff(y))

        fcoeffsrc.write(fsub_constr.format(
            name = routine_name,
            B = fprinter.doprint(C, 'C'),
            MB = fprinter.doprint(MC, 'MC'),
        ))

    #--------------------------------------------------------------------------#

    if nbl > 0:
        BL = Matrix(zeros((nbl,)))
        MBL = Matrix(zeros((nbl,ny)))

        for ib,b in zip(range(nbl),boundL):
            BL[ib] = discretize(-b)
            for iy,y in zip(range(ny),yvar):
                MBL[ib,iy] = discretize(b.diff(y))

        fcoeffsrc.write(fsub_bound.format(
            name = routine_name, lr = 'bl',
            B = fprinter.doprint(BL, 'B'),
            MB = fprinter.doprint(MBL, 'MB'),
        ))

    #--------------------------------------------------------------------------#

    if nbr > 0:
        BR = Matrix(zeros((nbr,)))
        MBR = Matrix(zeros((nbr,ny)))

        for ib,b in zip(range(nbr),boundR):
            BR[ib] = discretize(-b)
            for iy,y in zip(range(ny),yvar):
                MBR[ib,iy] = discretize(b.diff(y))

        fcoeffsrc.write(fsub_bound.format(
            name = routine_name, lr = 'br',
            B = fprinter.doprint(BR, 'B'),
            MB = fprinter.doprint(MBR, 'MB'),
        ))

    fcoeffsrc.write("\n")

    #--------------------------------------------------------------------------#

    for f in fswall: f.write('case({})\n'.format(model_nr))

    fswname.write ('  name = \"{}\"\n'.format(basename))

    fswptrs.write ('  fa  => {}\n'.format(routine_name))
    fswptrs.write ('  fc  => {}\n'      \
        .format(routine_name+'_c' if nc > 0 else 'NULL()'))
    fswptrs.write ('  fbl => {}\n'      \
        .format(routine_name+'_bl' if nbl > 0 else 'NULL()'))
    fswptrs.write ('  fbr => {}\n'      \
        .format(routine_name+'_br' if nbr > 0 else 'NULL()'))

    fswynum.write ('  ny = {}\n'.format(ny))

    fswdims.write ('  ny  = {}\n'.format(ny))
    fswdims.write ('  na  = {}\n'.format(na))
    fswdims.write ('  nc  = {}\n'.format(nc))
    fswdims.write ('  nbl = {}\n'.format(nbl))
    fswdims.write ('  nbr = {}\n'.format(nbr))

    fswhash.write ('  ihash = [{}]\n'.format(",".join([ str(i+1) \
            for i in yval_hash ])))

#------------------------------------------------------------------------------#

for f in fswall:
    f.write("case default\n  error stop\nend select\n")
    f.close()

fcoeffsrc.close()
