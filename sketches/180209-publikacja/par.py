# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from itertools import product as cartproduct

#------------------------------------------------------------------------------#

mbh, mdot, radius = 10.0, 0.05, 10.0
alpha, eta, nu = 0.1, 0.1, 0.5
# alpha, eta, nu = 0.1, 0.3, 0.0

#------------------------------------------------------------------------------#

figext = 'pdf'
# figext = 'png'

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

ix2fn = lambda s,ix: s[:-len(ix)] + ''.join([ '{}{:02d}'.format(a,i+1) for a,i in zip(s[-len(ix):],ix) ])
dict2cfg = lambda **kw: u''.join('{} {}\n'.format(k,v) for k,v in kw.items())

#------------------------------------------------------------------------------#

NN1 = 8

mdots_1 = np.logspace(-2.5, -0.5, NN1)
for i,mdot_ in enumerate(mdots_1):
    with open('data.1/01-' + ix2fn('M', (i,)) + '.par','w') as f:
        f.write(dict2cfg(mbh = mbh, mdot = mdot_, radius = radius,
            alpha = alpha, eta = eta, nu = nu))

alphas_1 = np.logspace(-2, 0, NN1) * eta
for i,alpha_ in enumerate(alphas_1):
    with open('data.1/01-' + ix2fn('A', (i,)) + '.par','w') as f:
        f.write(dict2cfg(mbh = mbh, mdot = mdot, radius = radius,
            alpha = alpha_, eta = eta, nu = nu))

nus_1 = np.logspace(-1, 1, NN1)
for i,nu_ in enumerate(alphas_1):
    with open('data.1/01-' + ix2fn('N', (i,)) + '.par','w') as f:
        f.write(dict2cfg(mbh = mbh, mdot = mdot, radius = radius,
            alpha = alpha, eta = eta, nu = nu_))

#------------------------------------------------------------------------------#

NN4 = 64
mdots = np.logspace(-2.0, -0.5, NN4)
alphas = np.linspace(1.0 / NN4, 1.0, NN4) * 0.3
radii = np.logspace(0.5, 2.0, NN4)
nus = np.logspace(-1, 1, NN4)

#------------------------------------------------------------------------------#

dsets = []

meta = {
    'M': (mdots, lbl_mdot, 'log'),
    'A': (alphas, lbl_alpha, 'linear'),
#    'Z': (zets, lbl_zet, 'linear'),
    'N': (nus, lbl_nu, 'log'),
    'R': (radii, lbl_radius, 'log'),
#    'B': (betas, lbl_beta, 'log'),
}

# MR MA MN NA

#------------------------------------------------------------------------------#
#
# ds = 'MR'
# dsets.append((ds, meta['M'], meta['R']))
# for i, mdot_ in enumerate(mdots):
#     for j, radius_ in enumerate(radii):
#         with open('data.0/' + ix2fn(ds,(i,j)) + '.par','w') as f:
#             f.write(dict2cfg(mbh = mbh, mdot = mdot_, radius = radius_,
#                 alpha = alpha, eta = eta, nu = nu))
#
# ds = 'MA'
# dsets.append((ds, meta['M'], meta['A']))
# for i, mdot_ in enumerate(mdots):
#     for j, alpha_ in enumerate(alphas):
#         with open('data.0/' + ix2fn(ds,(i,j)) + '.par','w') as f:
#             f.write(dict2cfg(mbh = mbh, mdot = mdot_, radius = radius,
#                 alpha = alpha_, eta = alpha_, nu = nu))
#
# ds = 'cMA'
# dsets.append((ds, meta['M'], meta['A']))
# for i, mdot_ in enumerate(mdots):
#     for j, alpha_ in enumerate(alphas):
#         with open('data.0/' + ix2fn(ds,(i,j)) + '.par','w') as f:
#             f.write(dict2cfg(mbh = mbh, mdot = mdot_, radius = radius,
#                 alpha = alpha_, eta = np.max(alphas), nu = nu))
#
# ds = 'NA'
# dsets.append((ds, meta['N'], meta['A']))
# for i, nu_ in enumerate(nus):
#     for j, alpha_ in enumerate(alphas):
#         with open('data.0/' + ix2fn(ds,(i,j)) + '.par','w') as f:
#             f.write(dict2cfg(mbh = mbh, mdot = mdot, radius = radius,
#                 alpha = alpha_, eta = alpha_, nu = nu_))
#
# ds = 'cNA'
# dsets.append((ds, meta['N'], meta['A']))
# for i, nu_ in enumerate(nus):
#     for j, alpha_ in enumerate(alphas):
#         with open('data.0/' + ix2fn(ds,(i,j)) + '.par','w') as f:
#             f.write(dict2cfg(mbh = mbh, mdot = mdot, radius = radius,
#                 alpha = alpha_, eta = np.max(alphas), nu = nu_))
#ds = 'NR'
#dsets.append((ds, meta['N'], meta['R']))
#for i, nu_ in enumerate(nus):
#    for j, radius_ in enumerate(radii):
#        with open('data.0/' + ix2fn(ds,(i,j)) + '.par','w') as f:
#            f.write(dict2cfg(mbh = mbh, mdot = mdot, radius = radius_,
#                alpha = alpha, eta = eta, nu = nu_))

#ds = 'AR'
#dsets.append((ds, meta['A'], meta['R']))
#for i, alpha_ in enumerate(alphas):
#    for j, radius_ in enumerate(radii):
#        with open('data.0/' + ix2fn(ds,(i,j)) + '.par','w') as f:
#            f.write(dict2cfg(mbh = mbh, mdot = mdot, radius = radius_,
#                alpha = alpha_, eta = alpha_, nu = nu))

ds = 'MN'
dsets.append((ds, meta['M'], meta['N']))
for i, mdot_ in enumerate(mdots):
    for j, nu_ in enumerate(nus):
        with open('data.0/' + ix2fn(ds,(i,j)) + '.par','w') as f:
            f.write(dict2cfg(mbh = mbh, mdot = mdot_, radius = radius,
                alpha = alpha, eta = eta, nu = nu_))

#------------------------------------------------------------------------------#

mdots_5 = np.round(np.logspace(np.log10(0.02), np.log10(0.06), 5), 4)
for i,mdot_ in enumerate(mdots_5):
    with open('data.1/05-' + ix2fn('M', (i,)) + '.par','w') as f:
        f.write(dict2cfg(mbh = mbh, mdot = mdot_, radius = radius,
            alpha = alpha, eta = eta, nu = nu))

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
