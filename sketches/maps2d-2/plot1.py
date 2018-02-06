from par import *

for sp in spaces:
    xx, xl, xs, xf = spaceinfo[sp[1]]
    yy, yl, ys, yf = spaceinfo[sp[0]]

    for sr in ['Q']:

        D = read2Ddataset(lambda i,j: fn(sp,i,j) + sr)

        fig, axes = plt.subplots(4, 5, sharex = True, sharey = True, figsize = (19,12.5))

        fig.suptitle(bigtitle)

        ax = axes[0,2]
        ax.set_title('coldens')
        zz = np.log10(DPA(D, 'coldens'))
        m = ax.pcolor(xx, yy, zz, cmap = 'YlOrBr')
        plt.colorbar(m, ax = ax)

        ax = axes[0,3]
        ax.set_title('instabil')
        zz = DPA(D, 'instabil')
        m = ax.pcolor(xx, yy, zz, cmap = 'PiYG', vmin = -0.1, vmax = 0.1)
        plt.colorbar(m, ax = ax)

        ax = axes[0,1]
        ax.set_title('xcor')
        zz = DPA(D, 'xcor')
        m = ax.pcolor(xx, yy, zz, cmap = 'RdGy', vmin = 0, vmax = 2)
        plt.colorbar(m, ax = ax)

        ax = axes[0,0]
        ax.set_title('beta_0')
        zz = DPA(D, 'beta_0')
        m = ax.pcolor(xx, yy, zz, cmap = 'BuPu_r', norm = LogNorm(vmin = 1e3, vmax = 1e0))
        plt.colorbar(m, ax = ax)

        ax = axes[0,4]
        ax.set_title('ddisk')
        zz = DPA(D, 'ddisk')
        m = ax.pcolor(xx, yy, zz, cmap = 'RdYlGn_r', vmin = 0, vmax = 0.02)
        plt.colorbar(m, ax = ax)

        for j, suff in zip(range(1,4), ['tmin', 'eqbc', 'therm']):

            ax = axes[j,0]
            ax.set_title('tavg_' + suff)
            zz = DPA(D, 'tavg_' + suff)
            m = ax.pcolor(xx, yy, zz, cmap = 'CMRmap', norm = LogNorm(vmin = 1e6, vmax = 1e8))
            plt.colorbar(m, ax = ax)

            ax = axes[j,1]
            ax.set_title('taues_' + suff)
            zz = DPA(D, 'taues_' + suff)
            m = ax.pcolor(xx, yy, zz, cmap = 'Spectral_r', norm = LogNorm(vmin = 3.16, vmax = 31.6)) # YlGnBu
            plt.colorbar(m, ax = ax)

            ax = axes[j,2]
            ax.set_title('mfrac_' + suff)
            zz = DPA(D, 'mfrac_' + suff)
            m = ax.pcolor(xx, yy, zz, cmap = 'YlOrBr', norm = LogNorm())
            plt.colorbar(m, ax = ax)

            ax = axes[j,3]
            ax.set_title('chi_' + suff)
            zz = DPA(D, 'chi_' + suff)
            m = ax.pcolor(xx, yy, zz, cmap = 'seismic', vmin = 0, vmax = 1)
            plt.colorbar(m, ax = ax)

            ax = axes[j,4]
            ax.set_title('compy_' + suff)
            zz = DPA(D, 'compy_' + suff) - DPA(D, 'compy_hard')
            m = ax.pcolor(xx, yy, zz, cmap = 'magma', vmin = 0)
            plt.colorbar(m, ax = ax)

        for axe in axes:
            for ax in axe:
                if sp == 'BN':
                    ax.plot(nus, alpha * (1 - nus),
                        color = '#ECD9EA', linestyle = '--')
                    ax.plot(nus, - alpha * (1 - nus),
                        color = '#ECD9EA', linestyle = ':')
                if sp == 'NB': ax.plot(alpha * (1 - nus), nus,
                        color = '#ECD9EA', linestyle = '--')
                if sp == 'BA': ax.plot(alphas, alphas * (1 - nu),
                        color = '#ECD9EA', linestyle = '--')
                if sp == 'AB': ax.plot(alphas * (1 - nu), alphas,
                        color = '#ECD9EA', linestyle = '--')
                if sp == 'BM': ax.axhline(alpha * (1 - nu),
                        color = '#ECD9EA', linestyle = '--')
                if sp == 'MB': ax.axvline(alpha * (1 - nu),
                        color = '#ECD9EA', linestyle = '--')
                if sp == 'AN':
                    ax.plot(1 + eta / alphas, alphas,
                        color = '#ECD9EA', linestyle = ':')
                    ax.plot(1 - eta / alphas, alphas,
                        color = '#ECD9EA', linestyle = '--')
                if sp == 'NA':
                    ax.plot(alphas, 1 + eta / alphas,
                        color = '#ECD9EA', linestyle = ':')
                    ax.plot(alphas, 1 - eta / alphas,
                        color = '#ECD9EA', linestyle = '--')
                ax.set_xlabel(xl)
                ax.set_xscale(xs)
                ax.set_xlim(np.min(xx), np.max(xx))
                ax.set_ylabel(yl)
                ax.set_yscale(ys)
                ax.set_ylim(np.min(yy), np.max(yy))

        plt.subplots_adjust(0.04, 0.05, 0.98, 0.93, 0.21, 0.25)
        plt.savefig('{sp}{sr}.png'.format(sp = sp, sr = sr))
        print '{sp}{sr}'.format(sp = sp, sr = sr)
        plt.close(fig)
