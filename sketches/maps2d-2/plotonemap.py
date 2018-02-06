from par import *
from sys import argv
from matplotlib.cm import coolwarm as cm

if len(argv) < 3:
    print u'usage: {} [field] [space1] [space2] ...'.format(argv[0])
    exit()

k = argv[1]

for sp in argv[2:]:

    xx, xl, xs, xf = spaceinfo[sp[1]]
    yy, yl, ys, yf = spaceinfo[sp[0]]
    sr = sp[2]

    zz = DPA(read2Ddataset(lambda i,j: fn(sp,i,j) + sr), k)

    fig, axes = plt.subplots(2, 2, figsize = (12,12))
    fig.suptitle(bigtitle)

    #---------------------------------------------------------------------------

    ax = axes[0,0]

    ax.set_xlabel(xl)
    ax.set_xscale(xs)
    ax.set_ylabel(yl)
    ax.set_yscale(ys)

    ax.set_title(k)
    m = ax.pcolor(xx, yy, zz, cmap = 'CMRmap')
    # plt.colorbar(m, ax = ax)

    #---------------------------------------------------------------------------

    ax = axes[1,0]

    ax.set_xlabel(xl)
    ax.set_xscale(xs)

    for i in range(0, NN, max(1,np.int(np.round(NN / 8.0)))):
        ax.plot(xx, zz[i,:], color = cm(i / (NN - 1.0)), label = u'{} = {}'.format(yl, yf.format(yy[i])))

    ax.legend(fontsize = 8)

    #---------------------------------------------------------------------------

    ax = axes[0,1]

    ax.set_ylabel(yl)
    ax.set_yscale(ys)

    for i in range(0, NN, max(1,np.int(np.round(NN / 8.0)))):
        ax.plot(zz[:,i], yy, color = cm(i / (NN - 1.0)), label = u'{:10} = {}'.format(xl, xf.format(xx[i])))

    ax.legend(fontsize = 8)


plt.show()
