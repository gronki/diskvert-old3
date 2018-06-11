# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from itertools import product as cartproduct

#------------------------------------------------------------------------------#

mbh, mdot, radius = 10.0, 0.05, 10
alpha, eta, nu = 0.1, 0.3, 0.5
# alpha, eta, nu = 0.1, 0.1, 1.0
# alpha, eta, nu = 0.1, 0.3, 0.0

#------------------------------------------------------------------------------#

# figext = 'pdf'
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
lbl_nu = u'$\\nu$' # u'$2 \\zeta / \\alpha_{{\\rm B}}$'
lbl_beta = u'$\\beta_0$'
lbl_zet = u'$\\xi$' # _{{\\rm rec}}

labels = {
    'M': lbl_mdot,
    'H': lbl_mbh,
    'R': lbl_radius,
    'A': lbl_alpha,
    'E': lbl_eta,
    'N': lbl_nu,
    'Z': lbl_zet,
    'B': lbl_beta,
}

#------------------------------------------------------------------------------#

ix2fn = lambda s,ix: s[:-len(ix)] + ''.join([ '{}{:02d}'.format(a,i+1) for a,i in zip(s[-len(ix):],ix) ])
dict2cfg = lambda **kw: u''.join('{:12s} {:.5g}\n'.format(k,v) for k,v in kw.items())
urange = lambda n: (np.arange(n, dtype = np.float64) + 1) / n

def mkpar(**kw):
    pars = dict(mbh = mbh, mdot = mdot, radius = radius, alpha = alpha,
        eta = eta, nu = nu)
    pars.update(kw)
    if 'alpha' in pars and 'zeta' in pars:
        pars['nu'] = 2 * pars['zeta'] / pars['alpha']
    return dict2cfg(**pars)

#------------------------------------------------------------------------------#

with open('data.2/cute.par','w') as f:
    f.write(mkpar(radius = 10, alpha = 0.05, eta = 0.5, zeta = 0.0))

#------------------------------------------------------------------------------#

NN1 = 7
dsets_multi = []
y = lambda x: (x + 2 * x**2) / 3.0

ds = 'M'
mdots_1 = np.logspace(-2.0, -0.3, NN1)
dsets_multi.append((ds, lbl_mdot, mdots_1))
for i, mdot_ in enumerate(mdots_1):
    with open('data.1/multi-' + ix2fn(ds, (i,)) + '.par', 'w') as f:
        f.write(mkpar(mdot = mdot_))

ds = 'R'
radii_1 = np.logspace(0.5, 2.0, NN1)
dsets_multi.append((ds, lbl_radius, radii_1))
for i, radius_ in enumerate(radii_1):
    with open('data.1/multi-' + ix2fn(ds, (i,)) + '.par', 'w') as f:
        f.write(mkpar(radius = radius_, mdot = 0.3))

alphas_1 = y(urange(NN1)) * eta

ds = 'etacon-nu0-A'
dsets_multi.append((ds, lbl_alpha, alphas_1))
for i, alpha_ in enumerate(alphas_1):
    with open('data.1/multi-' + ix2fn(ds, (i,)) + '.par', 'w') as f:
        f.write(mkpar(alpha = alpha_, eta = eta, nu = 0))

ds = 'etavar-nu0-A'
dsets_multi.append((ds, lbl_alpha, alphas_1))
for i, alpha_ in enumerate(alphas_1):
    with open('data.1/multi-' + ix2fn(ds, (i,)) + '.par', 'w') as f:
        f.write(mkpar(alpha = alpha_, eta = alpha_, nu = 0))

ds = 'etacon-nu1-A'
dsets_multi.append((ds, lbl_alpha, alphas_1))
for i, alpha_ in enumerate(alphas_1):
    with open('data.1/multi-' + ix2fn(ds, (i,)) + '.par', 'w') as f:
        f.write(mkpar(alpha = alpha_, eta = eta, nu = 1))

ds = 'etavar-nu1-A'
dsets_multi.append((ds, lbl_alpha, alphas_1))
for i, alpha_ in enumerate(alphas_1):
    with open('data.1/multi-' + ix2fn(ds, (i,)) + '.par', 'w') as f:
        f.write(mkpar(alpha = alpha_, eta = alpha_, nu = 1))

ds = 'N'
nus_1 = y(np.linspace(0, 1, NN1)) * (0.5 / alpha)
dsets_multi.append((ds, lbl_nu, nus_1))
for i,nu_ in enumerate(nus_1):
    with open('data.1/multi-' + ix2fn(ds, (i,)) + '.par', 'w') as f:
        f.write(mkpar(alpha = eta, eta = eta, nu = nu_))

#------------------------------------------------------------------------------#

NN4 = 24
dsets = []

#------------------------------------------------------------------------------#

mdots = np.logspace(-2.5, 0, NN4)
nus = np.logspace(-1, 1, NN4)

# ds = 'MN'
# dsets.append((ds, lbl_mdot, mdots, 'log', lbl_nu, nus, 'log'))
# for i, mdot_ in enumerate(mdots):
#     for j, nu_ in enumerate(nus):
#         with open('data.0/' + ix2fn(ds, (i,j)) + '.par', 'w') as f:
#             f.write(mkpar(mdot = mdot_, nu = nu_))

#------------------------------------------------------------------------------#

alphamax = 1.0

#------------------------------------------------------------------------------#

alphas = np.logspace(-2, np.log10(alphamax), NN4)

ds = 'etavar-nu0-MA'
dsets.append((ds, lbl_mdot, mdots, 'log', lbl_alpha, alphas, 'log'))
for i, mdot_ in enumerate(mdots):
    for j, alpha_ in enumerate(alphas):
        with open('data.0/' + ix2fn(ds, (i,j)) + '.par', 'w') as f:
            f.write(mkpar(mdot = mdot_, alpha = alpha_, eta = alpha_, nu = 0))

ds = 'etavar-nu1-MA'
dsets.append((ds, lbl_mdot, mdots, 'log', lbl_alpha, alphas, 'log'))
for i, mdot_ in enumerate(mdots):
    for j, alpha_ in enumerate(alphas):
        with open('data.0/' + ix2fn(ds, (i,j)) + '.par', 'w') as f:
            f.write(mkpar(mdot = mdot_, alpha = alpha_, eta = alpha_, nu = 1))

alphas = urange(NN4) * alphamax

ds = 'etacon-nu0-MA'
dsets.append((ds, lbl_mdot, mdots, 'log', lbl_alpha, alphas, 'linear'))
for i, mdot_ in enumerate(mdots):
    for j, alpha_ in enumerate(alphas):
        with open('data.0/' + ix2fn(ds, (i,j)) + '.par', 'w') as f:
            f.write(mkpar(mdot = mdot_, alpha = alpha_, eta = alphamax, nu = 0))

ds = 'etacon-nu1-MA'
dsets.append((ds, lbl_mdot, mdots, 'log', lbl_alpha, alphas, 'linear'))
for i, mdot_ in enumerate(mdots):
    for j, alpha_ in enumerate(alphas):
        with open('data.0/' + ix2fn(ds, (i,j)) + '.par', 'w') as f:
            f.write(mkpar(mdot = mdot_, alpha = alpha_, eta = alphamax, nu = 1))

#------------------------------------------------------------------------------#

radii = np.logspace(0.5, 2, NN4)

ds = 'MR'
dsets.append((ds, lbl_mdot, mdots, 'log', lbl_radius, radii, 'log'))
for i, mdot_ in enumerate(mdots):
    for j, radius_ in enumerate(radii):
        with open('data.0/' + ix2fn(ds, (i,j)) + '.par', 'w') as f:
            f.write(mkpar(mdot = mdot_, radius = radius_))

alphas = np.logspace(-2, np.log10(alphamax), NN4)

ds = 'etavar-nu0-AR'
dsets.append((ds, lbl_alpha, alphas, 'log', lbl_radius, radii, 'log'))
for i, alpha_ in enumerate(alphas):
    for j, radius_ in enumerate(radii):
        with open('data.0/' + ix2fn(ds, (i,j)) + '.par', 'w') as f:
            f.write(mkpar(radius = radius_, alpha = alpha_, eta = alpha_, nu = 0))

ds = 'etavar-nu1-AR'
dsets.append((ds, lbl_alpha, alphas, 'log', lbl_radius, radii, 'log'))
for i, alpha_ in enumerate(alphas):
    for j, radius_ in enumerate(radii):
        with open('data.0/' + ix2fn(ds, (i,j)) + '.par', 'w') as f:
            f.write(mkpar(radius = radius_, alpha = alpha_, eta = alpha_, nu = 1))

#------------------------------------------------------------------------------#

alphas = urange(NN4) * alphamax
zets = np.linspace(0, 0.5, NN4)

ds = 'etacon-AZ'
dsets.append((ds, lbl_alpha, alphas, 'linear', lbl_zet, zets, 'linear'))
for i, alpha_ in enumerate(alphas):
    for j, zet_ in enumerate(zets):
        with open('data.0/' + ix2fn(ds, (i,j)) + '.par', 'w') as f:
            f.write(mkpar(alpha = alpha_, eta = 0.5, zeta = zet_))

#------------------------------------------------------------------------------#

alphas = np.logspace(-2, np.log10(alphamax), NN4)

ds = 'etavar-AN'
dsets.append((ds, lbl_alpha, alphas, 'log', lbl_nu, nus, 'log'))
for i, alpha_ in enumerate(alphas):
    for j, nu_ in enumerate(nus):
        with open('data.0/' + ix2fn(ds, (i,j)) + '.par', 'w') as f:
            f.write(mkpar(alpha = alpha_, eta = alpha_, nu = nu_))

#------------------------------------------------------------------------------#

mdots_5 = np.array([0.025, 0.034, 0.045, 0.060])
for i,mdot_ in enumerate(mdots_5):
    with open('data.1/' + ix2fn('instabil-M', (i,)) + '.par', 'w') as f:
        f.write(mkpar(mbh = 10.0, mdot = mdot_, radius = 10.0,
            alpha = 0.1, eta = 0.1, nu = 1.0))

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
