from numpy import logspace, log10

fmt  = '  {:12s} {:10.5e}\n'
fmtc = '# {:12s} {:10.5e}\n'

mbh = 10.0
alpha = 0.007
betamax = (2 - alpha) / alpha

mdots = [0.0003, 0.0010, 0.0031, 0.0100, 0.0310]
betas = logspace(0, log10(betamax), 18)
rads =  logspace(log10(3.1), log10(160), 48)

izip = lambda x: zip(range(len(x)),x)

for imdot,mdot in izip(mdots):
    for ibeta,beta in izip(betas):
        zeta = alpha * (beta + 1) / 2
        for irad,rad in izip(rads):
            fn = 'M{:02d}R{:02d}B{:02d}'.format(imdot+1,irad+1,ibeta+1)
            with open('par/{}'.format(fn),'w') as f:
                f.write("# model: {}\n".format(fn))
                f.write( fmt.format('mbh',    mbh   ))
                f.write( fmt.format('mdot',   mdot  ))
                f.write( fmt.format('radius', rad   ))
                f.write( fmt.format('alpha',  alpha ))
                f.write( fmt.format('zeta',   zeta  ))
                f.write(fmtc.format('beta',   beta  ))

print "generated {} files".format(len(mdots)*len(betas)*len(rads))
