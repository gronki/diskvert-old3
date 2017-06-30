
from diskvert.cgs import *
from numpy import logspace, log10, sqrt
import matplotlib.pyplot as plt

trad, pgas, heat = 1.4837E+06, 1.4558E-02, 2.7746E+7

T = logspace(log10(trad),log10(trad)+3)

rho = pgas * 0.5 * cgs_mhydr / (cgs_boltz * T)

kappa_sct = 0.34
X0,Z0 = 0.72, 0.013
kabs0 = 3.68e22 * (1 - Z0) * (1 + X0) + 4.34e25 * Z0 * (1 + X0)
kappa_abs = lambda rho,T: kabs0 * rho * T**(-3.5)

cool = 4 * cgs_stef * rho * ( T - trad ) * (
            kappa_abs(rho,T) * (T + trad) * (T**2 + trad**2) \
            + kappa_sct * 4 * cgs_boltz * trad**4 / ( cgs_mel * cgs_c**2 )
        )

bil = cool/heat

plt.plot(rho,bil)
r = 1e-5
plt.ylim(-r,r)
plt.xscale('log')
plt.show()
