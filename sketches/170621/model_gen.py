import numpy as np
import numpy.linalg
from sympy import symbols, Symbol, Rational, Number, log as sympy_log, \
    expand_log, Matrix, exp as sympy_exp, nsimplify, solveset
from diskvert.cgs import *


mbh, mdot, r, alpha = symbols('mbh mdot r alpha', real = True, positive = True)
cgs_rschw_sun = 2 * cgs_graw * cgs_msun / cgs_c**2
cgs_mcrit_sun = 4 * cgs_pi * (cgs_graw * cgs_msun / cgs_c) \
* (cgs_mhydr / cgs_thomson)

fil = open('wzory_przyblizone.py','w')

logY_all = list()

for region in [1,2,3]:
    rho, T, H = symbols('rho{n} T{n} H{n}'.format(n = region), real = True, positive = True)

    r_cgs = r * mbh * cgs_rschw_sun
    mdot_cgs = mdot * mbh * cgs_mcrit_sun
    omega = ( cgs_graw * cgs_msun * mbh / r_cgs**3 )**0.5
    f = Symbol('f', real = True, positive = True)
    # kappa_ff = Symbol('kappa_ff', real = True, positive = True)
    kappa_ff = 6.13e22

    yvar = [rho, T, H]
    log_yvar = [ sympy_log(x) for x in yvar ]

    P1 = 4 * cgs_stef / (3 * cgs_c) * T**4
    P23 = 2 * cgs_boltz * rho * T / cgs_mhydr
    kap12 = cgs_kapes
    kap3 = kappa_ff * rho * T**Rational(-7,2)

    # in regions 1,2 electron scattering opacity, else free-free
    kap = kap3 if region == 3 else kap12
    # radiation pressure in the inner region, gas press everywhere else
    P =     P1 if region == 1 else P23

    # left and right hand sides of 3 equations
    L1 = 4 * cgs_pi * alpha * rho * H**3 * omega
    R1 = mdot_cgs * f
    L2 = ( Rational(3,4) * rho * H**2 )**2 * alpha * omega **3 * kap
    R2 = 4 * cgs_stef * T**4
    L3 = P
    R3 = rho * H**2 * omega**2

    # obtain so that RHS = 1 everywhere and we can logarithm
    A = [ L1/R1, L2/R2, L3/R3 ]
    # logarithm the lhs
    logA = [ expand_log(sympy_log(x)) for x in A ]

    # now collect the coefficients next to logarithms of each
    # variable: rho, T, H
    M = np.ndarray((3,3))
    for i in range(3):
        for j in range(3):
            M[i,j] = logA[i].coeff(log_yvar[j])
    print M

    # the amove will constitute our new left hand sides
    lhs = [ sum([ M[i,j] * log_yvar[j] for j in range(3) ]) for i in range(3) ]
    # right hand sides shall be a remainder
    rhs = [ lhs[i] - logA[i] for i in range(3) ]

    for i in range(3):
        print "{} = {}".format(lhs[i],rhs[i])

    # invert the matrix
    Minv = numpy.linalg.inv(M)
    # simplify the inversion coefficients before multiplication, e.g. 0.2 -> 1/5
    # calculate the variables
    logY = Matrix([ [ nsimplify(Minv[i,j]) for j in range(3) ] for i in range(3) ]) * Matrix(rhs)
    # add to the global array. we will need them later
    logY_all.append([logY[i] for i in range(3)])

    # for i in range(3):
    #     fil.write('log_{} = {}\n'.format(yvar[i], logY[i].evalf(5)))
    # fil.write('\n')
    for i in range(3):
        fil.write('{} = {}\n'.format(yvar[i], sympy_exp(logY[i]).evalf(5)))
    fil.write('\n')

temp = Symbol('logR')

fil.write("f12 = 1.0\nf23 = 1.0\n")

# przejscie 1->2
eq = nsimplify(logY_all[1][0]) - 3 * nsimplify(logY_all[1][1]) - sympy_log( (2 * cgs_stef * cgs_mhydr) / (3 * cgs_boltz * cgs_c) )
for x in solveset(eq.replace(sympy_log(r), temp),temp):
    # fil.write("log_r12 = {}\n".format(x.evalf(3)))
    fil.write("r12 = {}\n".format(sympy_exp(x.subs(f,'f12')).evalf(7)))

# przejscie 2->3
eq = nsimplify(logY_all[1][0]) - Rational(7,2) * nsimplify(logY_all[1][1]) - sympy_log( cgs_kapes / kappa_ff )
for x in solveset(eq.replace(sympy_log(r), temp),temp):
    # fil.write("log_r23 = {}\n".format(x.evalf(3)))
    fil.write("r23 = {}\n".format(sympy_exp(x.subs(f,'f23')).evalf(7)))

fil.write("f12 = 1 - sqrt(3/r12)\nf23 = 1 - sqrt(3/r23)\n")
