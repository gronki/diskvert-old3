import numpy as np
import matplotlib.pyplot as plt
from diskvert import *

clr = ['#93ABAA', '#09D6AB', '#CB2E13']

dd,pd = col2python('diskd.dat')
dw,pw = col2python('diskw.dat')
dc,pc = col2python('diskc.dat')

kabs = lambda rho, T: pd.kabp0 * rho * T**(-3.5)
yyc = pd.kappa_es * (4 * cgs_boltz * dd['trad']**4) / (cgs_mel * cgs_c**2)
tdumm1 = dd['trad'] + dd['heat'] / (4 * cgs_stef * dd['rho'] * yyc)
yyb = kabs(dd['rho'], tdumm1) * (dd['trad'] + tdumm1) * (dd['trad']**2 + tdumm1**2)
tdumm = dd['trad'] + dd['heat'] / (4 * cgs_stef * dd['rho'] * (yyb + yyc))
pdumm = dd['pgas'] * tdumm / dd['temp']

h0   = dc['h'][::2]
rho  = dc['rho'][::2]
heat = dc['heat'][::2]
trad = dc['trad'][::2]
temp = dc['temp'][::2]
T0 = np.logspace(np.log10(0.31 * pd.teff), np.log10(1e3 * pd.teff), 90)
hg, Tg = np.meshgrid(h0, T0)
bil = np.ndarray((len(T0),len(h0)))


for i in range(len(T0)):
    rho1 = rho
    # rho1 = rho * temp / T0[i]
    yyb = kabs(rho1, T0[i]) * (trad + T0[i]) * (trad**2 + T0[i]**2)
    yyc = pd.kappa_es * (4 * cgs_boltz * trad**4) / (cgs_mel * cgs_c**2)
    cool = 4 * cgs_stef * rho1 * (yyb + yyc) * (T0[i] - trad)
    bil[i,:] = 1 - cool / heat

plt.figure()
plt.pcolor(h0, T0, bil, cmap = 'bil1', vmin = -0.2, vmax = 0.2)
plt.loglog()
plt.xlim(2,2e2)

def plotthis(ax, col):
    ax.set_title(col)
    ax.plot(dd['h'], dd[col], color = clr[0], label = 'DYF')
    ax.plot(dw['h'], dw[col], color = clr[1], label = 'CMP')
    ax.plot(dc['h'], dc[col], color = clr[2], label = 'BIL')
    ax.set_xlim(2,2e2)
    ax.loglog()
    ax.legend()

fig, axes = plt.subplots(2, 2, figsize = (9,8), sharex = True)

plotthis(axes[0,0], 'rho')
plotthis(axes[0,1], 'pgas')
axes[0,1].plot(dd['h'], pdumm, '--', color = '#F450E3')
plotthis(axes[1,0], 'temp')
axes[1,0].plot(dd['h'], tdumm, '--', color = '#F450E3')
plotthis(axes[1,1], 'prad')

plt.show()
