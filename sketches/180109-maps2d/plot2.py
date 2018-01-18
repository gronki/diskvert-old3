from par import *
import matplotlib.pyplot as plt

for dset in ['W','C']:
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

    plt.subplot(3,4,1)
    plt.title(u'$\\tau_{\\rm es}$ of the corona')
    plt.pcolor(radii, etas, np.log10(DPA(D,'taues_cor')).transpose(), cmap = 'gist_earth_r')
    plt.loglog()
    plt.xlim(min(radii), max(radii))
    plt.xlabel('radius [x Rschw]')
    plt.ylabel('$\\eta = v_{\\rm B} / (\\Omega z)$')
    plt.colorbar()

    plt.subplot(3,4,5)
    plt.title(u'$\\tau_{\\rm es}$ of the therm. surface')
    plt.pcolor(radii, etas, np.log10(DPA(D,'taues_therm')).transpose(), cmap = 'gist_earth_r')
    plt.loglog()
    plt.xlim(min(radii), max(radii))
    plt.xlabel('radius [x Rschw]')
    plt.ylabel('$\\eta = v_{\\rm B} / (\\Omega z)$')
    plt.colorbar()

    plt.subplot(3,4,2)
    plt.title('$\\chi = F_{\\rm rad}^{\\rm cor} / F_{\\rm rad}^{\\rm tot}$')
    plt.pcolor(radii, etas, DPA(D,'chicor').transpose(), vmin = 0.2, vmax = 0.8, cmap = 'RdBu_r')
    plt.loglog()
    plt.xlim(min(radii), max(radii))
    plt.xlabel('radius [x Rschw]')
    plt.ylabel('$\\eta = v_{\\rm B} / (\\Omega z)$')
    plt.colorbar()

    plt.subplot(3,4,3)
    plt.title(u'$T_{\\rm av}$ to temp. minimum [K]')
    plt.pcolor(radii, etas, np.log10(DPA(D,'tcor')).transpose() / 11.6e6, cmap = 'CMRmap')
    plt.loglog()
    plt.xlim(min(radii), max(radii))
    plt.xlabel('radius [x Rschw]')
    plt.ylabel('$\\eta = v_{\\rm B} / (\\Omega z)$')
    plt.colorbar()

    plt.subplot(3,4,7)
    plt.title(u'column density $\\Sigma_{\\rm disk}$')
    plt.pcolor(radii, etas, (DPA(D,'coldens_disk').transpose()), cmap = 'CMRmap')
    plt.loglog()
    plt.xlim(min(radii), max(radii))
    plt.xlabel('radius [x Rschw]')
    plt.ylabel('$\\eta = v_{\\rm B} / (\\Omega z)$')
    plt.colorbar()

    plt.subplot(3,4,8)
    plt.title(u'column density $\\Sigma_{\\rm cor}$')
    # plt.title(u'tau thermal')
    plt.pcolor(radii, etas, DPA(D,'coldens_cor').transpose(), cmap = 'CMRmap')
    # plt.pcolor(radii, etas, DPA(D,'frfrac_cor').transpose(), cmap = 'RdBu_r', vmin = 0.2, vmax = 0.8)
    plt.loglog()
    plt.xlim(min(radii), max(radii))
    plt.xlabel('radius [x Rschw]')
    plt.ylabel('$\\eta = v_{\\rm B} / (\\Omega z)$')
    plt.colorbar()

    plt.subplot(3,4,4)
    plt.title(u'$T_{\\rm av}$ to thermaliz. depth [K]')
    plt.pcolor(radii, etas, np.log10(DPA(D,'tavg_therm')).transpose() / 11.6e6, cmap = 'CMRmap')
    plt.loglog()
    plt.xlim(min(radii), max(radii))
    plt.xlabel('radius [x Rschw]')
    plt.ylabel('$\\eta = v_{\\rm B} / (\\Omega z)$')
    plt.colorbar()


    plt.subplot(3,4,6)
    plt.title(u'$F_{\\rm mag}^{\\rm min} / F_{\\rm tot}^{\\rm min}$')
    plt.pcolor(radii, etas, DPA(D,'fbfrac_cor').transpose(), cmap = 'RdBu', vmin = 0.2, vmax = 0.8)
    plt.loglog()
    plt.xlim(min(radii), max(radii))
    plt.xlabel('radius [x Rschw]')
    plt.ylabel('$\\eta = v_{\\rm B} / (\\Omega z)$')
    plt.colorbar()

    #--------------------------------------------------------------

    plt.subplot(3,4,9)
    plt.title(u'tavg_eqbc')
    plt.yscale('log')
    plt.pcolor(radii, etas, np.log10(DPA(D, 'tavg_eqbc')).transpose(), cmap = 'CMRmap')
    plt.loglog()
    plt.xlim(min(radii), max(radii))
    plt.xlabel('radius [x Rschw]')
    plt.ylabel('$\\eta = v_{\\rm B} / (\\Omega z)$')
    plt.colorbar()

    plt.subplot(3,4,10)
    plt.title(u'taues_eqbc')
    plt.yscale('log')
    plt.pcolor(radii, etas, np.log10(DPA(D, 'taues_eqbc')).transpose(), cmap = 'gist_earth_r')
    plt.loglog()
    plt.xlim(min(radii), max(radii))
    plt.xlabel('radius [x Rschw]')
    plt.ylabel('$\\eta = v_{\\rm B} / (\\Omega z)$')
    plt.colorbar()

    plt.subplot(3,4,11)
    plt.title(u'compfr_tmin')
    plt.yscale('log')
    plt.pcolor(radii, etas, DPA(D, 'compfr_tmin').transpose(), cmap = 'CMRmap')
    plt.loglog()
    plt.xlim(min(radii), max(radii))
    plt.xlabel('radius [x Rschw]')
    plt.ylabel('$\\eta = v_{\\rm B} / (\\Omega z)$')
    plt.colorbar()

    plt.subplot(3,4,12)
    plt.title(u'compy_tmin')
    plt.yscale('log')
    plt.pcolor(radii, etas, DPA(D, 'compy_tmin').transpose(), cmap = 'CMRmap')
    plt.loglog()
    plt.xlim(min(radii), max(radii))
    plt.xlabel('radius [x Rschw]')
    plt.ylabel('$\\eta = v_{\\rm B} / (\\Omega z)$')
    plt.colorbar()

    #--------------------------------------------------------------

    plt.subplots_adjust(0.06, 0.09, 0.98, 0.90, 0.28, 0.32)

    fnpng = '{}-{}-2.png'.format(prefix,dset)
    plt.savefig(fnpng)
    print u'saved: {}'.format(fnpng)
