# coding: utf-8

from numpy import logspace,log10

mbh = 10
mdot = 0.007
radius = 9
alpha = 0.01

zeta = logspace(log10(alpha), 0, 48)
beta = 2 * zeta / alpha - 1

for bet,zet,i in zip(beta,zeta,range(len(zeta))):
    with open('B{:02d}.par'.format(i),'w') as f:
        f.write("mbh {}\nmdot {}\nradius {}\n".format(mbh,mdot,radius))
        f.write("alpha {}\nzeta {}\n".format(alpha,zet))
        f.write("# beta {}\n".format(bet))
