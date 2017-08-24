from numpy import logspace

alpha = 0.007
betamax = (2 - alpha) / alpha

mdots = logspace(-4,-1,7)
betas = logspace(0,np.log10(betamax),12)
rads =  logspace(0.5,2.7,72)

izip = lambda x: zip(range(len(x)),x)

for imdot,mdot in izip(mdots):
    for ibeta,beta in izip(betas):
        zeta = alpha * (beta + 1) / 2
        for irad,rad in izip(rads):
            fn = 'M{:02d}R{:02d}B{:02d}'.format(imdot+1,irad+1,ibeta+1)
            with open('{}.par'.format(fn),'w') as f:
                fmt = '  {:12s} {:10.5e}\n'
                fmtc = '# {:12s} {:10.5e}\n'
                f.write(fmt.format('mbh',10))
                f.write(fmt.format('mdot',mdot))
                f.write(fmt.format('radius',rad))
                f.write(fmt.format('alpha',alpha))
                f.write(fmt.format('zeta',zeta))
                f.write(fmtc.format('beta',beta))

print "generated {} files".format(len(mdots)*len(betas)*len(rads))
