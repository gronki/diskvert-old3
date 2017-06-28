# coding: utf-8
from fort import FortranProcedure, function_derivatives
from symbols import *
from sympy import Symbol, Function, Lambda, Derivative, IndexedBase, Eq, symbols

procedures = []

file_model_switch = open('src/model_select.inc','w')

# file_model_switch.write("submodule (relaxation) relaxation_model_select\n" \
#     "use relax_coefficients\nimplicit none\ncontains\n" \
file_model_switch.write("select case (model % nr)\n")

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

    # niewiadome w postaci listy
    yvar = [ rho, T_gas, F_rad ]
    if enableCorona: yvar.append(T_rad)
    if enableMagnetic: yvar.append(P_mag)
    if enableConduction: yvar.append(F_cond)

    yvar_all = [rho, T_gas, T_rad, F_rad, P_mag, F_cond]
    yval_hash = [ yvar.index(y) if y in yvar else -100 for y in yvar_all ]

    # nieprzezroczystosci i rpzewodnizctwa
    kappa_abs    = Function('fkabs')(rho,T_gas)
    kappa_cond   = Function('fkcnd')(rho,T_gas)
    global_functions = set([ kappa_abs, kappa_cond ])
    global_functions_extra = set([])

    fval = IndexedBase('fval', (3,len(global_functions)))

    # calkowita nieprzezroczystosc
    kappa = kappa_abs + kappa_es
    # predkosc wznoszenia pola
    vrise = Lambda(z, eta * Omega * z)
    # cisnienie gazu
    P_gas = k * rho * T_gas / ( mu * m_H )
    # cisnienie promieniowania
    P_rad = 4 * sigma / (3 * c) * T_rad**4
    # cisnienie calkowite
    P_tot = P_gas + P_rad + P_mag
    # strumien Poyntinga
    F_mag = 2 * P_mag * vrise(z)

    equations = []
    boundL = []
    boundR = []

    #
    # Hydrostatic equilibrium
    #
    equations.append(
        (P_gas.diff(z) - kappa * rho / c * F_rad + Derivative(P_mag,z))
        + Omega**2 * rho * z
    )

    #
    # Dyfuzja promieniowania
    #
    zeta_rad = 16 * sigma * T_rad**3 / (3 * kappa * rho)
    equations.append(
        F_rad + zeta_rad * Derivative(T_rad,z)
    )


    #
    # Bilans grzania i chlodzenia
    #

    GasHeating = alpha * Omega * P_tot
    if enableMagnetic: GasHeating = - vrise(z) * Derivative(P_mag,z)
    equations.append(
        Derivative(F_rad,z) + Derivative(F_cond,z) - GasHeating
    )
    boundL.append( F_rad )
    boundR.append( F_rad + F_mag + F_cond - F_acc )
    boundR.append( F_rad - 2 * sigma * T_rad**4 )

    #
    # Magnetic field evolution
    #
    if enableMagnetic:
        equations.append(
            2 * P_mag * vrise(z).diff(z) \
            + vrise(z) * Derivative(P_mag,z) \
            - alpha * Omega * P_tot
        )
        boundL.append( 2 * P_mag * vrise(z).diff(z) - alpha * Omega * P_tot )

    #
    # Funkcja zrodlowa
    #
    if enableCorona:
        SourceFunction = 4 * sigma * rho * ( T_gas - T_rad ) * (
            kappa_abs * (T_gas + T_rad) * (T_gas**2 + T_rad**2) \
            + kappa_es * 4 * k * T_rad**4 / ( m_el * c**2 )
        )
        equations.append(SourceFunction - Derivative(F_rad,z))
        boundL.append( T_gas - T_rad )

    #
    # Dyfuzja ciepla przez przewodnictwo
    #
    if enableConduction:
        equations.append(
            F_cond + kappa_cond * Derivative(T_gas,z)
        )
        boundL.append( F_cond )

    assert len(equations) == len(yvar)
    assert len(boundL) + len(boundR) == len(yvar)

    ny,neq,nbl,nbr = symbols('ny,neq,nbl,nbr', integer = True)
    neq = ny

    Y = IndexedBase('Y', (ny,))
    D = IndexedBase('D', (ny,))

    def discretize(eq):
        for iy,y in zip(range(len(yvar)),yvar):
            eq = eq.subs(Derivative(y,z),D[iy]).subs(y,Y[iy])
        eq, extrafun = function_derivatives(eq, global_functions)
        global_functions_extra.update(extrafun)
        return eq

    A = IndexedBase('A',(neq,))
    MY = IndexedBase('AY',(neq,ny))
    MD = IndexedBase('AD',(neq,ny))

    expr = []
    for ieq,eq in zip(range(len(equations)),equations):
        expr.append(Eq(A[ieq],discretize(eq)))
        for iy,y in zip(range(len(yvar)),yvar):
            expr.append(Eq(MY[ieq,iy],discretize(eq.diff(y))))
            expr.append(Eq(MD[ieq,iy],discretize(eq.diff(Derivative(y,z)))))

    expr_bl = []
    BL = IndexedBase('BL',(nbl,))
    MBL = IndexedBase('MBL',(nbl,ny))
    for ib,b in zip(range(len(boundL)),boundL):
        expr_bl.append(Eq(BL[ib],discretize(b)))
        for iy,y in zip(range(len(yvar)),yvar):
            expr_bl.append(Eq(MBL[ib,iy],discretize(b.diff(y))))

    expr_br = []
    BR = IndexedBase('BR',(nbr,))
    MBR = IndexedBase('MBR',(nbr,ny))
    for ib,b in zip(range(len(boundR)),boundR):
        expr_br.append(Eq(BR[ib],discretize(b)))
        for iy,y in zip(range(len(yvar)),yvar):
            expr_br.append(Eq(MBR[ib,iy],discretize(b.diff(y))))

    routine_name = 'COEFF_{magn}{comp}{cond}'.format(
        comp = 'COR' if enableCorona else 'DYF',
        magn = 'MAGN' if enableMagnetic else 'SS73',
        cond = 'CND' if enableConduction else '',
    )
    model_nr = 1 \
        + (1 if enableCorona else 0) \
        + (2 if enableMagnetic else 0) \
        + (4 if enableConduction else 0)
    file_model_switch.write('case({})\n'.format(model_nr))
    file_model_switch.write ('model % coeff => {}\n'.format(routine_name))
    file_model_switch.write ('model % ny = {}\n'.format(len(yvar)))
    file_model_switch.write ('model % nbl = {}\n'.format(len(boundL)))
    file_model_switch.write ('model % nbr = {}\n'.format(len(boundR)))
    file_model_switch.write ('model % coeff_L => {}\n'.format(routine_name+'_BL'))
    file_model_switch.write ('model % coeff_R => {}\n'.format(routine_name+'_BR'))
    file_model_switch.write ('model % ix = [{}]\n'.format(",".join([ str(i+1) for i in yval_hash ])))

    extern_functions = { f.func: str(f.func) for f in (global_functions | global_functions_extra) }

    procedures.append(FortranProcedure(
        routine_name, expr,
        kind = 'r64',
        arguments = [z, Y, D, A, MY, MD, ny],
        extern_functions = extern_functions,
        extern_variables = global_variables,
    ))

    # procedures.append(FortranProcedure(
    #     routine_name + '_SIZE', [
    #         Eq(ny,len(yvar)),
    #         Eq(nbl,len(boundL)),
    #         Eq(nbr,len(boundR)),
    #     ],
    #     arguments = [ny, nbl, nbr],
    # ))

    procedures.append(FortranProcedure(
        routine_name + '_BL', expr_bl,
        kind = 'r64',
        arguments = [z, Y, BL, MBL, ny, nbl],
        extern_functions = extern_functions,
        extern_variables = global_variables,
    ))

    procedures.append(FortranProcedure(
        routine_name + '_BR', expr_br,
        kind = 'r64',
        arguments = [z, Y, BR, MBR, ny, nbr],
        extern_functions = extern_functions,
        extern_variables = global_variables,
    ))

file_model_switch.write(
    "case default\n"
    "error stop\n"
    "end select\n"
)
# file_model_switch.write("end submodule\n")
file_model_switch.close()

with open('src/coefficients.f90','w') as f:
# from sys import stdout as f
    f.write("MODULE RELAX_COEFFICIENTS\n\n")
    f.write("USE IEEE_ARITHMETIC\n")
    f.write("USE SLF_CGS\n")
    f.write("use iso_fortran_env, only: r64 => real64\n")
    f.write("USE GLOBALS\n")
    f.write("\nIMPLICIT NONE\n")
    f.write("CONTAINS\n\n")
    for p in procedures:
        f.write(str(p) + "\n\n")
    f.write("\nEND MODULE\n")
