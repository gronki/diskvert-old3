import numpy as np
import matplotlib.pyplot as plt
from diskvert import *

mbh = 10
mdot = 0.06
radius = 7
radii = np.logspace(np.log10(3.5), np.log10(90), 9)
mdots = np.logspace(np.log10(0.03), np.log10(0.4), 9)
alpha = 0.01
eta = 0.5
etas = np.logspace(np.log10(alpha), 0, 9)
nu = 0
nus = np.logspace(-1, 2, 9)

#------------------------------------------------------------------------------#

par = """
mbh     {mbh}
mdot    {mdot}
alpha   {alpha}
radius  {radius}
zeta    {eta}
nu      {nu}
"""

#------------------------------------------------------------------------------#

if __name__ == '__main__':
    for i in range(len(etas)):
        with open('data/B{:02d}.par'.format(i+1), 'w') as f:
            f.write(par.format(
                mbh = mbh,
                mdot = mdot,
                alpha = alpha,
                radius = radius,
                eta = etas[i],
                nu = nu,
            ))
    for i in range(len(radii)):
        with open('data/R{:02d}.par'.format(i+1), 'w') as f:
            f.write(par.format(
                mbh = mbh,
                mdot = mdot,
                alpha = alpha,
                radius = radii[i],
                eta = eta,
                nu = nu,
            ))
    for i in range(len(mdots)):
        with open('data/M{:02d}.par'.format(i+1), 'w') as f:
            f.write(par.format(
                mbh = mbh,
                mdot = mdots[i],
                alpha = alpha,
                radius = radius,
                eta = eta,
                nu = nu,
            ))
    for i in range(len(nus)):
        with open('data/N{:02d}.par'.format(i+1), 'w') as f:
            f.write(par.format(
                mbh = mbh,
                mdot = mdot,
                alpha = alpha,
                radius = radius,
                eta = eta,
                nu = nus[i],
            ))
