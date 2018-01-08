from numpy import logspace, log10

fmt  = '  {:12s} {:10.5e}\n'
fmtc = '# {:12s} {:10.5e}\n'

mbh = 10.0
alpha = 0.01

mdots = [0.0003, 0.0010, 0.0031, 0.0100, 0.0310, 0.1000]
zetas = logspace(log10(alpha), 0, 16)
rads =  logspace(log10(3.1), log10(200), 48)

izip = lambda x: zip(range(len(x)), x)

fj = open('jobs.lst','w')

for imdot,mdot in izip(mdots):
    for izeta,zeta in izip(zetas):
        beta = 2 * zeta / alpha - 1
        for irad,rad in izip(rads):
            fn = 'M{:02d}R{:02d}B{:02d}'.format(imdot+1, irad+1, izeta+1)
            with open('par/{}'.format(fn),'w') as f:
                f.write("# model: {}\n".format(fn))
                f.write( fmt.format('mbh',    mbh   ))
                f.write( fmt.format('mdot',   mdot  ))
                f.write( fmt.format('radius', rad   ))
                f.write( fmt.format('alpha',  alpha ))
                f.write( fmt.format('zeta',   zeta  ))
                f.write(fmtc.format('beta',   beta  ))
            fj.write(fn + '\n')

fj.close()
print "generated {} files".format(len(mdots)*len(betas)*len(rads))
