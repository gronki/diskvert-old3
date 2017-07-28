from multiprocessing import cpu_count
import numpy as np

F = """#output {fn}
mbh {mbh:.5e}
mdot {mdot:.5e}
radius {radius:.5e}
alpha {alpha:.5e}
zeta {zeta:.5e}
# beta = {beta:.5e}
"""

nmdot = 10
nbeta = 7
nrad = 18

izip = lambda x: zip(range(len(x)),x)
alpha = 0.01
for imdot,mdot in izip(np.logspace(-3,-1,nmdot)):
    betamax = 2 / alpha - 1
    for ibeta,beta in izip(np.logspace(0,np.log10(betamax),nbeta)):
        zeta = alpha * (beta + 1) / 2
        for irad,rad in izip(np.logspace(np.log10(5),np.log10(50),nrad)):
            fn = 'M{:02d}R{:02d}B{:02d}'.format(imdot,irad,ibeta)
            with open('{}.par'.format(fn),'w') as f:
                f.write(F.format(
                    fn = fn,
                    mbh = 10.0,
                    mdot = mdot,
                    radius = rad,
                    alpha = alpha,
                    zeta = zeta,
                    beta = beta,
                ))
print "generated {} files".format(nmdot*nbeta*nrad)
