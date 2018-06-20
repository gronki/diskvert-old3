from random import uniform

pars = dict(mbh = 10, mdot = uniform(0.045, 0.055), radius = uniform(5.0, 7.0),
    alpha = uniform(0.015, 0.025), eta = 0.20)

for k,v in pars.iteritems():
    print "{:10s} {:.5g}".format(k,v)
