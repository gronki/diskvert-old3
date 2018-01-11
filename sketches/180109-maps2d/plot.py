from par import *
from diskvert import pyminiconf

dataset = 'dataw'
D = read2Ddataset(nrad, N2, lambda i,j: dataset + '/' + fn(i,j) + '.txt')

import matplotlib.pyplot as plt

plt.figure(figsize = (12,8))

plt.subplot(231)
plt.title(u'$\\tau_{\\rm es}$ of the corona')
plt.pcolor(radii, etas, DPA(D,'taues_cor').transpose(), cmap = 'CMRmap_r')
plt.loglog()
plt.xlim(min(radii), max(radii))
plt.xlabel('radius [x Rschw]')
plt.ylabel('$v_{\\rm mag} / (\\Omega z)$')
plt.colorbar()

plt.subplot(232)
plt.title('energy in the corona')
plt.pcolor(radii, etas, DPA(D,'chicor').transpose(), vmin = 0.2, vmax = 0.8, cmap = 'RdBu_r')
plt.loglog()
plt.xlim(min(radii), max(radii))
plt.xlabel('radius [x Rschw]')
plt.ylabel('$v_{\\rm mag} / (\\Omega z)$')
plt.colorbar()

plt.subplot(233)
plt.title(u'$T_{\\rm av}$ of the corona [keV]')
plt.pcolor(radii, etas, DPA(D,'tcor_keV').transpose(), cmap = 'CMRmap')
plt.loglog()
plt.xlim(min(radii), max(radii))
plt.xlabel('radius [x Rschw]')
plt.ylabel('$v_{\\rm mag} / (\\Omega z)$')
plt.colorbar()

from matplotlib.cm import coolwarm as cm

plt.subplot(234)
dd = DPA(D, 'taues_cor', mask = False)
for i in range(0,N2):
    plt.plot(radii, dd[:,i], color = cm((i + 1.0) / N2), label = '$\\eta = {:.3f}, \\beta = {:.2e}$'.format(etas[i], betas[i]))
plt.xscale('log')
plt.xlim(min(radii), max(radii))
plt.xlabel('radius [x Rschw]')
plt.ylabel('$\\tau_{\\rm es}$ of the corona')

plt.subplot(235)
dd = DPA(D, 'chicor', mask = False)
for i in range(0,N2):
    plt.plot(radii, dd[:,i], color = cm((i + 1.0) / N2), label = '$\\eta = {:.3f}, \\beta = {:.2e}$'.format(etas[i], betas[i]))
plt.xscale('log')
plt.xlim(min(radii), max(radii))
plt.xlabel('radius [x Rschw]')
plt.ylabel('energy in the corona')

plt.subplot(236)
dd = DPA(D, 'tcor_keV', mask = False)
for i in range(0,N2):
    plt.plot(radii, dd[:,i], color = cm((i + 1.0) / N2), label = '$\\eta = {:.3f}, \\beta = {:.2e}$'.format(etas[i], betas[i]))
plt.xscale('log')
plt.xlim(min(radii), max(radii))
plt.xlabel('radius [x Rschw]')
plt.ylabel('$T_{\\rm av}$ of the corona [keV]')

# plt.tight_layout()
plt.subplots_adjust(0.06, 0.07, 0.98, 0.95, 0.29, 0.21)
plt.show()
