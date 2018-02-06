import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from diskvert import *


#------------------------------------------------------------------------------#

mbh     = 10.0
mdot    = 0.05
radius  = 5.50
alpha   = 0.10
eta     = 0.50
nu      = 0.00

#------------------------------------------------------------------------------#

NN = 30

#------------------------------------------------------------------------------#

mbhs   = np.logspace(np.log10(1.0), np.log10(1e9), NN)
mdots  = np.logspace(np.log10(0.003), np.log10(0.5), NN)
radii  = np.logspace(np.log10(3.8), np.log10(50), NN)

alphas = np.logspace(np.log10(eta) - 2, np.log10(eta), NN)
# etas   = alpha * (np.logspace(-2, np.log10(2 / alpha - 1), NN) + 1) / 2
etas = np.logspace(np.log10(0.55 * alpha), 0, NN)
nus    = np.logspace(np.log10(0.1), np.log10(10.0), NN)

#------------------------------------------------------------------------------#

fn = lambda s,i,j: "{A}{i:02d}{B}{j:02d}".format(A = s[0], B = s[1], i = i + 1, j = j + 1)

#------------------------------------------------------------------------------#

def DPA(D, k, mask = True):
    from numpy import vectorize
    okay = lambda d: d != None and hasattr(d,k) and type(getattr(d,k)) != str and d.converged
    d0 = vectorize(lambda d: getattr(d,k) if okay(d) else np.nan)(D)
    if not mask: return d0
    from numpy.ma import masked_where
    dm = vectorize(lambda d: not okay(d))(D)
    return masked_where(dm, d0)

#------------------------------------------------------------------------------#

def read2Ddataset(fn):
    from diskvert import pyminiconf
    from numpy import ndarray
    D = ndarray((NN,NN), dtype = object)
    for i in range(NN):
        for j in range(NN):
            try:
                with open('data/' + fn(i,j) + '.txt', 'r') as f:
                    D[i,j] = pyminiconf(f)
            except IOError:
                D[i,j] = None
    return D

#------------------------------------------------------------------------------#

spaces = [ 'NH', 'NM', 'NR', 'NA' ]
spaces = [ 'MH', 'NH', 'BH', 'BM', 'NM', 'BR', 'MR', 'AR', 'HR', 'NR', 'BN' ]
spaces = [ 'AM', 'BM', 'NM', 'AN', 'BN', 'AB' ]
spaces = [ 'BM', 'BN' ]

spaceinfo = {
    'M': (mdots, '$\\dot{{m}}$', 'log', '{:.4f}'),
    'B': (etas, '$\\eta$', 'log', '{:.3f}'),
    'H': (mbhs, '$M_{{\\rm BH}}$', 'log', '{:.2e}'),
    'N': (nus, '$\\nu_{{\\rm rec}}$', 'log', '{:.1f}'),
    'A': (alphas, '$\\alpha_{{\\rm B}}$', 'log', '{:.3f}'),
    'R': (radii, '$R / R_{{\\rm schw}}$', 'log', '{:.1f}'),
}

bigtitle = u'$\\log M_{{\\rm BH}}$ = {:.1f}, $\\dot{{m}}$ = {:.5f}, $R / R_{{\\rm schw}}$ = {:.1f}, $\\alpha_{{\\rm B}}$ = {:.3f}, $\\eta$ = {:.3f}, $\\nu_{{\\rm rec}}$ = {:.1f}'.format(np.log10(mbh), mdot, radius, alpha, eta, nu)

#------------------------------------------------------------------------------#

if __name__ == '__main__':

    par =  u"mbh    {mbh:.4e}\n"    \
            "mdot   {mdot:.6f}\n"   \
            "alpha  {alpha:.6f}\n"  \
            "radius {radius:.3f}\n" \
            "zeta   {eta:.6f}\n"    \
            "nu     {nu:.3f}\n"

    fnlist = []

    def gen(s):
        fnl = []
        for i in range(NN):
            for j in range(NN):
                with open('data/' + fn(s,i,j) + '.par','w') as f:
                    f.write(par.format(
                    mbh =    mbhs[i]    if s[0] == 'H'  \
                        else mbhs[j]    if s[1] == 'H'  \
                        else mbh,
                    mdot =   mdots[i]   if s[0] == 'M'  \
                        else mdots[j]   if s[1] == 'M'  \
                        else mdot,
                    alpha =  alphas[i]  if s[0] == 'A'  \
                        else alphas[j]  if s[1] == 'A'  \
                        else alpha,
                    eta =    etas[i]    if s[0] == 'B'  \
                        else etas[j]    if s[1] == 'B'  \
                        else eta,
                    nu =     nus[i]     if s[0] == 'N'  \
                        else nus[j]     if s[1] == 'N'  \
                        else nu,
                    radius = radii[i]   if s[0] == 'R'  \
                        else radii[j]   if s[1] == 'R'  \
                        else radius,
                    ))
                    fnl.append(fn(s,i,j))
        return fnl

    for sp in spaces: fnlist.extend(gen(sp))

    with open('data/jobs.lst','w') as f:
        for x in fnlist: f.write(x + '\n')
