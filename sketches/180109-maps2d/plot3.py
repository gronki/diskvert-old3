from par import *
import matplotlib.pyplot as plt

for dset in ['W','Wp','C']:
    D = read2Ddataset(nrad, N2, lambda i,j: fn(i,j) + dset + '.txt')

    plt.figure(figsize = (16,11))
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

    from matplotlib.cm import coolwarm_r as cm

    def plotser(ax, k):
        dd = DPA(D, k, mask = False)
        for i in range(0, nrad, 3):
            ax.plot(etas, dd[i,:], color = cm((i + 1.0) / nrad), label = u'$R = {:.1f}$'.format(radii[i]))
        ax.set_xscale('log')
        ax.set_xlim(min(etas), max(etas))
        ax.set_xlabel('$\\eta = v_{\\rm B} / (\\Omega z)$')

    plt.subplot(3,4,1)
    plt.title(u'$\\tau_{\\rm es}$ at $T_{{\\rm min}}')
    plotser(plt.gca(), 'taues_tmin')
    plt.yscale('log')
    # plt.ylim(2e2,5)
    plt.ylabel('$\\tau_{\\rm es}$')

    plt.subplot(3,4,5)
    plt.title(u'$\\tau_{\\rm es}$ at $\\tau* = 1$')
    plotser(plt.gca(), 'taues_therm')
    plt.yscale('log')
    # plt.ylim(2e2,5)
    plt.ylabel('$\\tau_{\\rm es}$')

    plt.subplot(3,4,2)
    plt.title('$\\chi = F_{\\rm rad}^{\\rm cor} / F_{\\rm rad}^{\\rm tot}$')
    plotser(plt.gca(), 'chi_tmin')
    plt.ylim(0,1)
    plt.ylabel('$\\chi$')

    plt.subplot(3,4,3)
    plt.title(u'$T_{\\rm av}$ at $T_{{\\rm min}}$ [K]')
    plt.yscale('log')
    plt.ylabel('$T_{\\rm av}$')
    # plt.ylim(temp_estim, temp_estim * 30)
    plotser(plt.gca(), 'tavg_tmin')


    plt.subplot(3,4,4)
    plt.title(u'$T_{\\rm av}$ at $\\tau* = 1$ [K]')
    plt.yscale('log')
    plt.ylabel('$T_{\\rm av}$')
    # plt.ylim(temp_estim, temp_estim * 30)
    plotser(plt.gca(), 'tavg_therm')


    plt.subplot(3,4,7)
    plt.title(u'column density $\\Sigma_{\\rm disk}$')
    plt.yscale('log')
    plt.ylabel('$\\Sigma_{\\rm disk}$')
    plotser(plt.gca(), 'coldens_below_tmin')

    plt.subplot(3,4,8)
    plt.title(u'column density $\\Sigma_{\\rm cor}$')
    plt.yscale('log')
    plt.ylabel('$\\Sigma_{\\rm cor}$')
    plotser(plt.gca(), 'coldens_tmin')

    plt.subplot(3,4,6)
    plt.title(u'$F_{\\rm mag} / F_{\\rm tot}$ at $T_{{\\rm min}}$')
    plt.ylim(0,1)
    plotser(plt.gca(), 'fbfrac_tmin')

    #--------------------------------------------------------------

    plt.subplot(3,4,9)
    plt.title(u'tavg_eqbc')
    plt.yscale('log')
    plotser(plt.gca(), 'tavg_eqbc')

    plt.subplot(3,4,10)
    plt.title(u'taues_eqbc')
    plt.yscale('log')
    plotser(plt.gca(), 'taues_eqbc')

    plt.subplot(3,4,11)
    plt.title(u'compfr_tmin')
    plt.yscale('log')
    plotser(plt.gca(), 'compfr_tmin')

    plt.subplot(3,4,12)
    plt.title(u'compy_tmin')
    plt.yscale('log')
    plotser(plt.gca(), 'compy_tmin')

    #--------------------------------------------------------------

    plt.subplots_adjust(0.06, 0.09, 0.98, 0.90, 0.28, 0.32)

    fnpng = '{}-{}-3.png'.format(prefix,dset)
    plt.savefig(fnpng)
    print u'saved: {}'.format(fnpng)
