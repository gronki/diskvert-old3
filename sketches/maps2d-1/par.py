import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from diskvert import *

kv0 = [1, 4, 3, 1]

from sys import argv
if len(argv) == 5:
    kv = [ int(x) if x not in '?X' else kv0[i] for i,x in zip(range(4),argv[1:5]) ]
else: kv = kv0
#------------------------------------------------------------------------------#

mbh     = [10.0, 1e7, 1e8][kv[0] - 1]
mdot    = [0.0003, 0.0010, 0.0030, 0.0100, 0.0300, 0.1000, 0.3000][kv[1] - 1]
alpha   = [0.003, 0.010, 0.030, 0.100, 0.300][kv[2] - 1]
nu      = [0.0, 1.0, 3.0, 9.0][kv[3] - 1]

print u"mbh = {:.2e}, mdot = {}, alpha = {}, nu = {}".format(mbh, mdot, alpha, nu)
#------------------------------------------------------------------------------#

nrad = 40
N2 = nrad

#------------------------------------------------------------------------------#

radii = np.logspace(np.log10(4.0), np.log10(200.0), nrad)
etas = np.logspace(np.log10(alpha), 0, N2)
betas = 2 * etas / alpha + nu - 1

#------------------------------------------------------------------------------#

series = 'AGN' if mbh > 1e3 else 'XRB'
prefix = '{}-{:04.1f}-M{:06.4f}-A{:05.3f}-N{:04.1f}'.format(series, np.log10(mbh), mdot, alpha, nu)
fn = lambda i,j: "R{:02d}B{:02d}".format(i+1, j+1)

#------------------------------------------------------------------------------#

def DPA(D, k, mask = True):
    from numpy import vectorize
    okay = lambda d: d != None and hasattr(d,k) and d.converged
    d0 = vectorize(lambda d: getattr(d,k) if okay(d) else np.nan)(D)
    if not mask: return d0
    from numpy.ma import masked_where
    dm = vectorize(lambda d: not okay(d))(D)
    return masked_where(dm, d0)

#------------------------------------------------------------------------------#

def read2Ddataset(nx, ny, fnpath):
    from diskvert import pyminiconf
    from numpy import ndarray
    D = ndarray((nx,ny), dtype = object)
    for i in range(nx):
        for j in range(ny):
            try:
                with open('data.' + prefix + '/' + fnpath(i,j), 'r') as f:
                    D[i,j] = pyminiconf(f)
            except IOError:
                D[i,j] = None
    return D

#------------------------------------------------------------------------------#

if __name__ == '__main__':
    par = """
    mbh     {mbh}
    mdot    {mdot}
    alpha   {alpha}
    radius  {radius}
    zeta    {eta}
    nu      {nu}
    """

    from os import mkdir
    mkdir('data.' + prefix)
    print prefix

    with open('data.{}/run.sh'.format(prefix),'w') as f:
        f.write(u'cat jobs.lst | parallel --bar bash ../job.sh\n')

    fnlist = []

    for i in range(nrad):
        for j in range(N2):
            with open('data.' + prefix + '/' + fn(i,j) + '.par','w') as f:
                f.write(par.format(
                    mbh = mbh,
                    mdot = mdot,
                    alpha = alpha,
                    radius = radii[i],
                    eta = etas[j],
                    nu = nu,
                ))
                fnlist.append(fn(i,j))

    with open('data.{}/jobs.lst'.format(prefix),'w') as f:
        for x in fnlist: f.write(x + '\n')
