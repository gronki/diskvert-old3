from par import *
from matplotlib.cm import coolwarm as cm
from matplotlib.cm import coolwarm_r as cm_r
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------#

model = 'C'

#------------------------------------------------------------------------------#

def prepaxes():
    fig, axes = plt.subplots(2, 3, figsize = (14,9))

    axes[0,0].loglog()
    axes[0,0].set_xlim(2,2e2)
    axes[0,0].set_xlabel('z / H')
    axes[0,0].set_title('T [K]')
    axes[0,0].set_ylim(7e5,4e9)

    axes[0,1].set_xscale('log')
    axes[0,1].set_xlim(2,2e2)
    axes[0,1].set_xlabel('z / H')
    axes[0,1].set_yscale('log')
    axes[0,1].set_ylim(1e-8,3)
    axes[0,1].set_title('density [cgs]')

    axes[0,2].set_xscale('log')
    axes[0,2].set_xlim(2,2e2)
    axes[0,2].set_xlabel('z / H')
    axes[0,2].set_yscale('log')
    axes[0,2].set_ylim(1e3,1e-4)
    axes[0,2].set_title('magnetic $\\beta = P_{gas} / P_{mag}$')

    axes[1,1].set_xscale('log')
    axes[1,1].set_xlim(2,2e2)
    axes[1,1].set_xlabel('z / H')
    axes[1,1].set_yscale('log')
    axes[1,1].set_ylim(1e4,1e-3)
    axes[1,1].set_title('optical depth $\\tau_{{\\rm es}}$')

    axes[1,2].set_xscale('log')
    axes[1,2].set_xlim(1e3,1e-3)
    axes[1,2].set_xlabel('$\\tau_{{\\rm es}}$')
    axes[1,2].set_title('Compton cooling fraction $\\Lambda_C / (\\Lambda_B + \\Lambda_C)$')

    axes[1,0].set_xscale('log')
    axes[1,0].set_yscale('log')
    axes[1,0].set_xlim(1e3,1e-3)
    axes[1,0].set_xlabel('$\\tau_{{\\rm es}}$')
    axes[1,0].set_ylim(7e5,4e9)
    axes[1,0].set_title('T [K]')

    return fig, axes

#------------------------------------------------------------------------------#

fig, axes = prepaxes()

N = len(etas)
for i in range(N):
    d,p = col2python('data/B{:02d}{}.dat'.format(i+1, model))
    lbl = '$\\eta$ = {:.2f}'.format(etas[i])
    if not p.converged: continue

    axes[0,0].plot(d['h'], d['temp'], color = cm((i + 1.0) / N), label = lbl)
    axes[1,0].plot(d['taues'], d['temp'], color = cm((i + 1.0) / N), label = lbl)
    axes[0,1].plot(d['h'], d['rho'], color = cm((i + 1.0) / N), label = lbl)
    axes[0,2].plot(d['h'], d['beta'], color = cm((i + 1.0) / N), label = lbl)
    axes[1,1].plot(d['h'], d['taues'], color = cm((i + 1.0) / N), label = lbl)
    axes[1,2].plot(d['taues'], d['coolc'] / (d['coolb'] + d['coolc']), color = cm((i + 1.0) / N), label = lbl)
    axes[1,2].axhline(0.5, color = '#348339', linewidth = 0.7)


axes[0,0].legend(fontsize = 9, loc = 'best')
plt.tight_layout()
plt.savefig('pl1.png')


#------------------------------------------------------------------------------#

fig, axes = prepaxes()

N = len(radii)
for i in range(N):
    d,p = col2python('data/R{:02d}{}.dat'.format(i+1, model))
    lbl = '$R$ = {:.2f}'.format(radii[i])
    if not p.converged: continue


    axes[0,0].plot(d['h'], d['temp'], color = cm_r((i + 1.0) / N), label = lbl)
    axes[1,0].plot(d['taues'], d['temp'], color = cm_r((i + 1.0) / N), label = lbl)
    axes[0,1].plot(d['h'], d['rho'], color = cm_r((i + 1.0) / N), label = lbl)
    # coolr = d['coolc'] / (d['coolb'] + d['coolc'])
    axes[0,2].plot(d['h'], d['beta'], color = cm_r((i + 1.0) / N), label = lbl)
    axes[1,1].plot(d['h'], d['taues'], color = cm_r((i + 1.0) / N), label = lbl)
    axes[1,2].plot(d['taues'], d['coolc'] / (d['coolb'] + d['coolc']), color = cm_r((i + 1.0) / N), label = lbl)
    axes[1,2].axhline(0.5, color = '#348339', linewidth = 0.7)


axes[0,0].legend(fontsize = 9, loc = 'best')

plt.tight_layout()
plt.savefig('pl2.png')

#------------------------------------------------------------------------------#

fig, axes = prepaxes()

N = len(mdots)
for i in range(N):
    d,p = col2python('data/M{:02d}{}.dat'.format(i+1, model))
    lbl = '$\\dot m$ = {:.3f}'.format(mdots[i])
    if not p.converged: continue

    axes[0,0].plot(d['h'], d['temp'], color = cm((i + 1.0) / N), label = lbl)
    axes[1,0].plot(d['taues'], d['temp'], color = cm((i + 1.0) / N), label = lbl)
    axes[0,1].plot(d['h'], d['rho'], color = cm((i + 1.0) / N), label = lbl)
    # coolr = d['coolc'] / (d['coolb'] + d['coolc'])
    axes[0,2].plot(d['h'], d['beta'], color = cm((i + 1.0) / N), label = lbl)
    axes[1,1].plot(d['h'], d['taues'], color = cm((i + 1.0) / N), label = lbl)
    axes[1,2].plot(d['taues'], d['coolc'] / (d['coolb'] + d['coolc']), color = cm((i + 1.0) / N), label = lbl)
    axes[1,2].axhline(0.5, color = '#348339', linewidth = 0.7)


axes[0,0].legend(fontsize = 9, loc = 'best')
plt.tight_layout()
plt.savefig('pl3.png')

#------------------------------------------------------------------------------#

fig, axes = prepaxes()

N = len(nus)
for i in range(N):
    d,p = col2python('data/N{:02d}{}.dat'.format(i+1, model))
    lbl = '$\\nu$ = {:.1f}'.format(nus[i])
    if not p.converged: continue

    axes[0,0].plot(d['h'], d['temp'], color = cm((i + 1.0) / N), label = lbl)
    axes[1,0].plot(d['taues'], d['temp'], color = cm((i + 1.0) / N), label = lbl)
    axes[0,1].plot(d['h'], d['rho'], color = cm((i + 1.0) / N), label = lbl)
    axes[0,2].plot(d['h'], d['beta'], color = cm((i + 1.0) / N), label = lbl)
    axes[1,1].plot(d['h'], d['taues'], color = cm((i + 1.0) / N), label = lbl)
    axes[1,2].plot(d['taues'], d['coolc'] / (d['coolb'] + d['coolc']), color = cm((i + 1.0) / N), label = lbl)
    axes[1,2].axhline(0.5, color = '#348339', linewidth = 0.7)


axes[0,0].legend(fontsize = 9, loc = 'best')
plt.tight_layout()
plt.savefig('pl4.png')

#------------------------------------------------------------------------------#
