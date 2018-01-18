from par import *
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

for dset in ['W','C']:
    D = read2Ddataset(nrad, N2, lambda i,j: fn(i,j) + dset + '.txt')

    def plotser(ax, k):
        from matplotlib.cm import coolwarm as cm
        dd = DPA(D, k, mask = False)
        for i in range(0, N2, 3):
            ax.plot(radii, dd[:,i], color = cm((i + 1.0) / nrad), label = u'$R = {:.1f}$'.format(radii[i]))
        ax.set_xscale('log')
        ax.set_xlim(min(radii), max(radii))
        ax.set_xlabel('$R / R_{\\rm schw}$')

    plt.figure(figsize = (16,8))
    plt.suptitle('$\\log M_{{\\rm BH}} = {logM:.1f}$, $\\dot m = {mdot:.4f}$, $\\alpha_B = {alpha:.4f}$, $\\nu = {nu:.1f}$, $R = [{rin:.1f},{rout:.1f}] \\times R_{{\\rm schw}}$'.format(
        logM = np.log10(mbh),
        mdot = mdot,
        alpha = alpha,
        nu = nu,
        rin = min(radii),
        rout = max(radii),
    ))

    temp_estim = 0.26213 * (cgs_c**5 * mdot / (cgs_graw * cgs_kapes \
        * cgs_stef * mbh * sol_mass))**0.25

    plt.subplot(2,4,1)
    plt.title(u'$\\tau_{\\rm es}$ at $T_{{\\rm min}}$')
    plotser(plt.gca(), 'taues_cor')
    plt.yscale('log')
    plt.ylim(0.1,100)

    plt.subplot(2,4,3)
    plt.title('$\\chi = F_{\\rm rad}^{\\rm cor} / F_{\\rm rad}^{\\rm tot}$')
    plotser(plt.gca(), 'chicor')
    plt.ylim(0.2,0.8)

    plt.subplot(2,4,2)
    plt.title(u'$T_{\\rm av}$ at $T_{{\\rm min}}$')
    plt.ylim(3e5,5e8)
    plotser(plt.gca(), 'tcor')
    plt.yscale('log')

    plt.subplot(2,4,4)
    plt.title(u'$T_{\\rm eff}$')
    plt.ylim(3e5,5e8)
    plotser(plt.gca(), 'teff')
    plt.yscale('log')

    plt.subplot(2,4,7)
    plt.title(u'column density $\\Sigma_{\\rm disk}$')
    plotser(plt.gca(), 'coldens_disk')
    plt.yscale('log')

    plt.subplot(2,4,8)
    plt.title(u'column density $\\Sigma_{\\rm cor}$')
    plotser(plt.gca(), 'coldens_cor')
    plt.yscale('log')


    #--------------------------------------------------------------

    plt.subplot(2,4,6)
    plt.title(u'$T_{\\rm avg}$ at $\\Lambda_{{\\rm C}} = \\Lambda_{{\\rm B}}$')
    plotser(plt.gca(), 'tavg_eqbc')
    plt.ylim(3e5,5e8)
    plt.yscale('log')

    plt.subplot(2,4,5)
    plt.title(u'$\\tau_{\\rm es}$ at $\\Lambda_{{\\rm C}} = \\Lambda_{{\\rm B}}$')
    plotser(plt.gca(), 'taues_eqbc')
    plt.ylim(0.1,100)
    plt.yscale('log')


    #--------------------------------------------------------------

    plt.subplots_adjust(0.06, 0.09, 0.98, 0.90, 0.30, 0.32)

    fnpng = '{}-{}-6b.png'.format(prefix,dset)
    plt.savefig(fnpng)
    print u'saved: {}'.format(fnpng)
