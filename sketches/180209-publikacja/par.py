# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from itertools import product as cartproduct

#------------------------------------------------------------------------------#

mbh, mdot, radius = 10.0, 0.05, 6.0
alpha, eta, nu = 0.1, 0.1, 1.0

#------------------------------------------------------------------------------#

figext = 'pdf'
figext = 'png'

figsc = 0.65
figw_sm = figsc * 8.8
figw_md = figsc * 12.0
figw_lg = figsc * 17.0

#------------------------------------------------------------------------------#

lbl_mbh = u'$M_{{\\rm BH}}$'
lbl_mdot = u'$\\dot{{m}}$'
lbl_radius = u'$R / R_{{\\rm schw}}$'
lbl_alpha = u'$\\alpha_{{\\rm B}}$'
lbl_eta = u'$\\eta$'
lbl_nu = u'$\\nu_{{\\rm rec}}$'
lbl_zet = u'$\\zeta_{{\\rm rec}}$'
lbl_beta = u'$\\beta_0$'

#------------------------------------------------------------------------------#

ix2fn = lambda s,ix: ''.join([ '{}{:02d}'.format(a,i+1) for a,i in zip(s,ix) ])
dict2cfg = lambda **kw: u''.join('{} {}\n'.format(k,v) for k,v in kw.items())

#------------------------------------------------------------------------------#

NN1 = 8

mdots_1 = np.logspace(-2.5, -0.5, NN1)
for i,x in enumerate(mdots_1):
    with open('data.1/01-' + ix2fn('M', (i,)) + '.par','w') as f:
        f.write(dict2cfg(mbh = mbh, mdot = x, radius = radius,
            alpha = alpha, zeta = eta, nu = nu))

alphas_1 = np.logspace(-2, 0, NN1) * eta
for i,x in enumerate(alphas_1):
    with open('data.1/01-' + ix2fn('A', (i,)) + '.par','w') as f:
        f.write(dict2cfg(mbh = mbh, mdot = mdot, radius = radius,
            alpha = x, zeta = eta, nu = nu))

nus_1 = np.logspace(-1, 2, NN1)
for i,x in enumerate(alphas_1):
    with open('data.1/01-' + ix2fn('N', (i,)) + '.par','w') as f:
        f.write(dict2cfg(mbh = mbh, mdot = mdot, radius = radius,
            alpha = alpha, zeta = eta, nu = x))

#------------------------------------------------------------------------------#

NN4 = 40
mdots = np.logspace(-2.5, -0.5, NN4)
alphas = np.linspace(0.4 / NN4, 0.4, NN4)
radii = np.logspace(0.5, 2, NN4)
zets = np.linspace(0, 0.7, NN4)
nus = np.logspace(-1, 1, NN4)
betas = np.logspace(0, 1.5, NN4)

dsets = []

#------------------------------------------------------------------------------#

dsets.append('MR')
for i,j in cartproduct(range(NN4), range(NN4)):
    with open('data.0/04-' + ix2fn('MR',(i,j)) + '.par','w') as f:
        f.write(dict2cfg(mbh = mbh, mdot = mdots[i], radius = radii[j],
            alpha = alpha, zeta = eta, nu = nu))

#------------------------------------------------------------------------------#

dsets.append('MA')
for i,j in cartproduct(range(NN4), range(NN4)):
    with open('data.0/04-' + ix2fn('MA',(i,j)) + '.par','w') as f:
        f.write(dict2cfg(mbh = mbh, mdot = mdots[i], radius = radius,
            alpha = alphas[j], zeta = alphas[j], nu = nu))

#------------------------------------------------------------------------------#

dsets.append('MN')
for i,j in cartproduct(range(NN4), range(NN4)):
    with open('data.0/04-' + ix2fn('MN',(i,j)) + '.par','w') as f:
        f.write(dict2cfg(mbh = mbh, mdot = mdots[i], radius = radius,
            alpha = alpha, zeta = eta, nu = nus[j]))

#------------------------------------------------------------------------------#

dsets.append('NA')
for i,j in cartproduct(range(NN4), range(NN4)):
    with open('data.0/04-' + ix2fn('NA',(i,j)) + '.par','w') as f:
        f.write(dict2cfg(mbh = mbh, mdot = mdot, radius = radius,
            alpha = alphas[j], zeta = alphas[j], nu = nus[i]))

#------------------------------------------------------------------------------#

meta = {
    'M': (mdots, lbl_mdot, 'log'),
    'A': (alphas, lbl_alpha, 'linear'),
    'Z': (zets, lbl_zet, 'linear'),
    'N': (nus, lbl_nu, 'log'),
    'R': (radii, lbl_radius, 'log'),
    'B': (betas, lbl_beta, 'log'),
}

#------------------------------------------------------------------------------#

mdots_5 = np.round(np.logspace(np.log10(0.018), np.log10(0.040), 5), 4)
for i,x in enumerate(mdots_5):
    with open('data.1/05-' + ix2fn('M', (i,)) + '.par','w') as f:
        f.write(dict2cfg(mbh = mbh, mdot = x, radius = radius,
            alpha = alpha, zeta = eta, nu = nu))

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
