# coding: utf-8

from numpy import logspace,log10

mbh = 10
mdot = 0.007
radius = 9
alpha = 0.01

beta = logspace(log10(0.333), log10(2 - alpha) - log10(alpha), 36)
zeta = alpha * (beta + 1) / 2

fn = lambda a,b: 'B{:06d}'.format(int(round(b*100)))
for bet,zet in zip(beta,zeta):
    with open(fn(alpha,bet) + '.par','w') as f:
        f.write("mbh {}\nmdot {}\nradius {}\n".format(mbh,mdot,radius))
        f.write("alpha {}\nzeta {}\n".format(alpha,zet))
        f.write("# beta {}\n".format(bet))
