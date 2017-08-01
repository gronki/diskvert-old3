from multiprocessing import cpu_count
import numpy as np

nmdot = 5
nbeta = 18
nrad = 36
alpha = 0.01
alpha = 0.007

izip = lambda x: zip(range(len(x)),x)

for imdot,mdot in izip(np.logspace(-3,-1,nmdot)):
    betamax = (2 - alpha) / alpha
    for ibeta,beta in izip(np.logspace(0,np.log10(betamax),nbeta)):
        zeta = alpha * (beta + 1) / 2
        for irad,rad in izip(np.logspace(np.log10(3.1),np.log10(50),nrad)):
            fn = 'M{:03d}R{:03d}B{:03d}'.format(imdot+1,irad+1,ibeta+1)
            with open('{}.par'.format(fn),'w') as f:
                fmt = '  {:12s} {:10.5e}\n'
                fmtc = '# {:12s} {:10.5e}\n'
                f.write(fmt.format('mbh',10))
                f.write(fmt.format('mdot',mdot))
                f.write(fmt.format('radius',rad))
                f.write(fmt.format('alpha',alpha))
                f.write(fmt.format('zeta',zeta))
                f.write(fmtc.format('beta',beta))

print "generated {} files".format(nmdot*nbeta*nrad)
