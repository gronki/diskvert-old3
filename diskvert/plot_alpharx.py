# coding: utf-8

import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib import rc


def read_alpharx(n, name = 'disk'):
    from numpy import loadtxt
    return loadtxt('{0}.{1:03}.dat'.format(name,n),
        dtype = dict(
            names = ('n','z','rho','T','Trad','Frad','beta','Fcond'),
            formats = ('i4','f4','f4','f4','f4','f4','f4','f4')
        ),
    )

def axprep(plt):
    fig, axes = plt.subplots(1, 3, figsize = (12,4))
    ax = axes[0]
    ax.set_title('rho')
    ax.set_yscale('log')
    ax = axes[1]
    ax.set_title('T')
    ax = axes[2]
    ax.set_title('Frad')
    return fig, axes

def axplot(axes, data, **kwargs):
    axes[0].plot(data['z'], data['rho'], **kwargs)
    axes[1].plot(data['z'], data['Trad'], linewidth = 0.8, **kwargs)
    axes[1].plot(data['z'], data['T'], **kwargs)
    axes[2].plot(data['z'], data['Frad'], **kwargs)

def parseargv():
    from sys import argv
    name = 'disk'
    frames = [ int(arg) for arg in argv[1:] ]
    return name,frames

def main_plotmulti():
    name, frames = parseargv()
    for frame in frames:
        data = read_alpharx(frame, name)
        fig, axes = axprep(plt)
        axplot(axes,data)
        plt.savefig('{0}.{1:03}.png'.format(name,frame), dpi = 144)
        plt.close(fig)
    print 'plotted {} frames'.format(len(frames))

def main_plotcumul():
    name, frames = parseargv()
    cmap = get_cmap('inferno_r')
    fig, axes = axprep(plt)
    for frame in frames:
        data = read_alpharx(frame, name)
        c = cmap(1.0 * frame / max(frames))
        axplot(axes, data, color = c)
    plt.savefig('{0}multi{1}.png'.format(name,len(frames)), dpi = 144)
    plt.close(fig)
    print 'combined {} plots'.format(len(frames))

if __name__ == '__main__':
    main_plotmulti()
