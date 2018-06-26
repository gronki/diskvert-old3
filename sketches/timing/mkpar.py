from random import uniform

pars = dict(mbh = uniform(9.5, 10.5), mdot = uniform(0.045, 0.060), radius = uniform(4.0, 8.0),
    alpha = uniform(0.015, 0.025), eta = 0.20)

for k,v in pars.iteritems():
    print "{:10s} {:.5g}".format(k,v)
