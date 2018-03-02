from diskvert import *
import numpy as np
import matplotlib.pyplot as plt
from par import *

beta = np.logspace(-0.5,2.5)

plt.figure()

# plt.plot(beta, 0.25 * beta**-0.55)
plt.plot(beta, 2 * 0.1 / (1 + beta), color = '#F9CEB9', label = '$\\eta = 0.1$')
plt.plot(beta, 2 * 0.3 / (1 + beta), color = '#E07A5E', label = '$\\eta = 0.3$')
plt.plot(beta, 2 * 0.9 / (1 + beta), color = '#D70000', label = '$\\eta = 0.9$')
# plt.plot(beta, 2 * 0.1 / (1 + beta - 0.77), '--', color = '#B3D5F5')
# plt.plot(beta, 2 * 0.3 / (1 + beta - 0.77), '--', color = '#68AFF1')
# plt.plot(beta, 2 * 0.9 / (1 + beta - 0.77), '--', color = '#0070D7')
plt.loglog()
plt.legend()

from symresults import *

plt.scatter(sym_betas, sym_alphas, color = '#D70000')
plt.errorbar(sym_betas, sym_alphas, xerr = sym_betas_err, yerr = sym_alphas_err, color = '#D70000', fmt = 'o')

# sym_betas_2 = np.array([68, 41, 17, 2.1, 0.4, 0.31])
# sym_alphas_2 = np.array([0.0088, 0.014, 0.035, 0.3, 1.0, 1.1])
# plt.scatter(sym_betas_2, sym_alphas_2, color = '#0070D7')

plt.xlabel('$\\beta_0$')
plt.ylabel('$\\alpha_{{\\rm B}}$')


# sym_etas = sym_alphas * (sym_betas + 1) / 2
# print 'etas -> ', sym_etas
# print ' eta_m = ', np.mean(sym_etas), np.std(sym_etas)
#
# sym_betas = np.array([68, 41, 17, 2.1, 0.4, 0.31])
# sym_alphas = np.array([0.0088, 0.014, 0.035, 0.3, 1.0, 1.1])
# plt.figure()
# nus = np.linspace(-1,2,250)
# err = [np.std(sym_alphas * (sym_betas + 1 - nu) / 2) for nu in nus]
# plt.plot(nus,err)
#
# sym_etas = sym_alphas * (sym_betas + 1 - 0.77) / 2
# print 'etas -> ', sym_etas
# print ' eta_m = ', np.mean(sym_etas), np.std(sym_etas)

plt.savefig('3.{}'.format(figext))

exit()

plt.figure()
plt.ylabel('$\\eta$')
plt.xlabel('$\\alpha_{{\\rm B}}$')
plt.scatter(sym_alphas, sym_etas)
a = sum(sym_etas * sym_alphas) / sum(sym_alphas**2)
x = np.linspace(0,0.4)
plt.plot(x, a*x, '--')
print 'eta = {} * alpha'.format(a)
plt.savefig('3b.png')

plt.figure()
plt.scatter(sym_betas, (2 * (sym_etas + sym_zetas) / sym_alphas - 1) / sym_betas - 1)
plt.savefig('3c.png')
