from par import *
import matplotlib.pyplot as plt

for dset in ['W','C']:
    D = read2Ddataset(nrad, N2, lambda i,j: fn(i,j) + dset + '.txt')

    plt.figure(figsize = (12,8))

    plt.subplot(231)
    plt.title(u'$\\tau_{\\rm es}$ of the corona')
    plt.pcolor(radii, etas, DPA(D,'taues_cor').transpose(), cmap = 'CMRmap_r')
    plt.loglog()
    plt.xlim(min(radii), max(radii))
    plt.xlabel('radius [x Rschw]')
    plt.ylabel('$\\eta = v_{\\rm B} / (\\Omega z)$')
    plt.colorbar()

    plt.subplot(232)
    plt.title('$\\chi = F_{\\rm rad}^{\\rm cor} / F_{\\rm rad}^{\\rm tot}$')
    plt.pcolor(radii, etas, DPA(D,'chicor').transpose(), vmin = 0.2, vmax = 0.8, cmap = 'RdBu_r')
    plt.loglog()
    plt.xlim(min(radii), max(radii))
    plt.xlabel('radius [x Rschw]')
    plt.ylabel('$\\eta = v_{\\rm B} / (\\Omega z)$')
    plt.colorbar()

    plt.subplot(233)
    plt.title(u'$T_{\\rm av}$ in the corona [keV]')
    plt.pcolor(radii, etas, DPA(D,'tcor_keV').transpose(), cmap = 'CMRmap')
    plt.loglog()
    plt.xlim(min(radii), max(radii))
    plt.xlabel('radius [x Rschw]')
    plt.ylabel('$\\eta = v_{\\rm B} / (\\Omega z)$')
    plt.colorbar()

    from matplotlib.cm import coolwarm as cm

    skip = 3

    plt.subplot(234)
    dd = DPA(D, 'taues_cor', mask = False)
    for i in range(0,N2,skip):
        plt.plot(radii, dd[:,i], color = cm((i + 1.0) / N2), label = '$\\eta = {:.3f}, \\beta = {:.2e}$'.format(etas[i], betas[i]))
    plt.xscale('log')
    plt.xlim(min(radii), max(radii))
    plt.xlabel('radius [x Rschw]')
    plt.ylabel('$\\tau_{\\rm es}$ of the corona')

    plt.subplot(235)
    dd = DPA(D, 'chicor', mask = False)
    for i in range(0,N2,skip):
        plt.plot(radii, dd[:,i], color = cm((i + 1.0) / N2), label = '$\\eta = {:.3f}, \\beta = {:.2e}$'.format(etas[i], betas[i]))
    plt.xscale('log')
    plt.xlim(min(radii), max(radii))
    plt.xlabel('radius [x Rschw]')
    plt.ylabel('$\\chi = F_{\\rm rad}^{\\rm cor} / F_{\\rm rad}^{\\rm tot}$')

    plt.subplot(236)
    dd = DPA(D, 'tcor_keV', mask = False)
    for i in range(0,N2,skip):
        plt.plot(radii, dd[:,i], color = cm((i + 1.0) / N2), label = '$\\eta = {:.3f}, \\beta = {:.2e}$'.format(etas[i], betas[i]))
    plt.xscale('log')
    plt.xlim(min(radii), max(radii))
    plt.xlabel('radius [x Rschw]')
    plt.ylabel('$T_{\\rm av}$ in the corona [keV]')

    plt.subplots_adjust(0.06, 0.07, 0.98, 0.95, 0.29, 0.21)

    fnpng = '{}-{}.png'.format(prefix,dset)
    plt.savefig(fnpng)
    print u'saved: {}'.format(fnpng)
