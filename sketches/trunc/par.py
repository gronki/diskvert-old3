import numpy as np
import matplotlib.pyplot as plt
from itertools import product as cartproduct

#------------------------------------------------------------------------------#

mbh = 10
mdot = 0.05
radius = 10
alpha = 0.1
eta = 0.3
nu = 0.0

#------------------------------------------------------------------------------#

lbl_mbh = u'$M_{{\\rm BH}}$'
lbl_mdot = u'$\\dot{{m}}$'
lbl_radius = u'$R / R_{{\\rm schw}}$'
lbl_alpha = u'$\\alpha_{{\\rm B}}$'
lbl_eta = u'$\\eta$'
lbl_nu = u'$\\nu_{{\\rm rec}}$'

#------------------------------------------------------------------------------#

ix2fn = lambda s,ix: ''.join([ '{}{:03d}'.format(a,i+1) for a,i in zip(s,ix) ])
dict2cfg = lambda **kw: u''.join('{} {}\n'.format(k,v) for k,v in kw.items())

#------------------------------------------------------------------------------#

nus = (1 + 18)**np.linspace(0, 1, 120) - 1
zetas = np.linspace(0,1,120)

mdots = np.logspace(-2.5, -0.5, 6)
for i,j in cartproduct(range(len(mdots)), range(len(nus))):
    with open('data/' + ix2fn('MN',(i,j)) + '.par','w') as f:
        f.write(dict2cfg(mbh = mbh, mdot = mdots[i], radius = radius,
            alpha = alpha, zeta = eta, nu = 2 * zetas[j] / alpha))

alphas = np.logspace(-2, 0, 6) * eta
for i,j in cartproduct(range(len(alphas)), range(len(nus))):
    with open('data/' + ix2fn('AN',(i,j)) + '.par','w') as f:
        f.write(dict2cfg(mbh = mbh, mdot = mdot, radius = radius,
            alpha = alphas[i], zeta = eta, nu = 2 * zetas[j] / alphas[i]))

# for i,j in cartproduct(range(len(alphas)), range(len(nus))):
#     with open('data/' + ix2fn('AN',(i,j)) + '.par','w') as f:
#         f.write(dict2cfg(mbh = mbh, mdot = mdot, radius = radius,
#             alpha = alphas[i], zeta = 0.7 * alphas[i], nu = nus[j]))

#------------------------------------------------------------------------------#

def read2Ddataset(fn, N1, N2):
    from diskvert import pyminiconf
    from numpy import ndarray
    dset = ndarray((N1,N2), dtype = object)
    for i in range(N1):
        for j in range(N2):
            try:
                with open(fn(i,j), 'r') as f:
                    dset[i,j] = pyminiconf(f)
            except IOError:
                dset[i,j] = None
    return dset

#------------------------------------------------------------------------------#

def DPA(dset, k, mask = True):

    okay = lambda d: d != None \
        and hasattr(d,k) and type(getattr(d,k)) != str \
        and hasattr(d,'converged') and d.converged

    from numpy import vectorize
    d0 = vectorize(lambda d: getattr(d,k) if okay(d) else np.nan)(dset)
    if not mask: return d0

    from numpy.ma import masked_where
    dm = vectorize(lambda d: not okay(d))(dset)
    return masked_where(dm, d0)
