from par import *
from matplotlib.cm import coolwarm as cm
from matplotlib.cm import coolwarm_r as cm_r
import matplotlib.pyplot as plt

fig,axes = plt.subplots(1, 2, figsize = (10,4.5))

for i in range(5, N2, 6):
    d,p = col2python('data.' + prefix + '/' + fn(20,i) + 'W.dat')
    lbl = '$\\eta$ = {:.2f}'.format(etas[i])

    axes[0].plot(d['h'], d['temp'], color = cm((i + 1.0) / N2), label = lbl)

    coolr = d['coolc'] / (d['coolb'] + d['coolc'])
    axes[1].plot(d['taues'], coolr, color = cm((i + 1.0) / N2), label = lbl)

axes[0].loglog()
axes[0].set_xlim(2,2e2)
axes[0].legend(fontsize = 9, loc = 'best')
axes[0].set_xlabel('z / H')
axes[0].set_ylabel('T [K]')

axes[1].set_xscale('log')
axes[1].set_xlim(1e3,1e-3)
axes[1].set_xlabel('z / H')
axes[1].set_ylabel('T [K]')

plt.show()
