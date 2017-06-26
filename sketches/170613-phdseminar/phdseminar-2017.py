#!/usr/bin/env python
# coding: utf-8

from diskvert import *
import numpy as np

# initialize exapmle accretion disk
dv_init_disk(10,0.003,10)
dv_init_ss73(0.05)
dv_eval_globals()

zscal = 6.4651E+04

# prepare memtory
nz = 2**10
z  = np.ndarray(nz)
y  = np.ndarray((4,nz),  order = 'F')
dy = np.ndarray((4,nz),  order = 'F')
a  = np.ndarray((8,nz),  order = 'F')

# run integration
dv_run_ss73(z,nz,y,dy,a)

# interpolate initial values for new model
nx = 300
ny = 3
x = np.linspace(0, z.max(), nx)
Y = np.ndarray(nx*ny)
# Y[0::ny] = np.interp(x, z[::-1], a[2,::-1])
# Y[1::ny] = np.interp(x, z[::-1], a[1,::-1])
# Y[2::ny] = np.interp(x, z[::-1], y[2,::-1])
Y[0::ny] = np.linspace(1e-3, 0, nx)
Y[1::ny] = np.linspace(2.71 * a[1,0], a[1,0], nx)
Y[2::ny] = np.linspace(0, y[2,0], nx)

T_floor = a[1,0]

# matrix
A = np.zeros((nx*ny))
M = np.zeros((nx*ny,nx*ny), order = 'F')

# change parameters a little bit
# dv_init_ss73(0.05)
# dv_eval_globals()

# draw
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from matplotlib.cm import get_cmap
from matplotlib import rc

cmap = get_cmap('YlOrRd')
cmap = get_cmap('Greys')



def ramp(it, niter, n = 1, ramp_0 = 0):
    t = np.clip(it / (niter - 1.0), 0, 1)
    return ramp_0 + (1 - ramp_0) * t**n

def ramp_smooth(it, niter, ramp0 = 0):
    t = np.clip(it / (niter - 1.0), 0, 1)
    return ramp0 + (1 - ramp0) * (3 * t**2 - 2 * t**3)

def ramp_trig(it, niter, ramp0 = 0):
    t = np.clip(it / (niter - 1.0), 0, 1)
    return ramp0 + (1 - ramp0) * (1 - np.cos(np.pi * t)) / 2

niter = 128
niter_full = 8

rc('font', family='serif', size=11)

dY = np.ndarray(nx*ny)
special = [0, niter-1]
for it in range(niter):
    fig,ax = plt.subplots(1, 3, figsize=(16 * 0.7,9 * 0.7))

    plt.subplots_adjust(left=0.05, right=0.98, top=0.94, bottom=0.08)

    color_gray = '#999999'

    ax[0].set_title('$\\rho$ [${\\rm g} / {\\rm cm}^3$]')
    ax[0].set_xlabel('z / H')
    ax[0].set_ylim([0,5e-2])
    ax[0].plot(z / zscal, a[2,:], color=color_gray)

    ax[1].set_title('$T$ [K]')
    ax[1].set_xlabel('z / H')
    ax[1].set_ylim([0,1e7])
    ax[1].plot(z / zscal, a[1,:], color=color_gray)

    ax[2].set_title('$F_{\\rm rad}$')
    ax[2].set_xlabel('z / H')
    ax[2].set_ylim([0,2e20])
    ax[2].plot(z / zscal, y[2,:], color=color_gray)

    # generate the matrix and solve
    dv_generate_corrections(x,nx,Y,dY,ny)

    rmp = ramp_trig(it, niter - niter_full, 0.01)

    err = np.sqrt(np.sum( np.nan_to_num(dY)**2 ))
    print "{:5d} {:8.2e} {:5.3f}".format(it,err,rmp)

    c = 'black' #if it in special else cmap(ramp(it, niter, 0.5))
    lthck = 1.25 #if it in special else 0.8

    ax[0].plot(x / zscal, Y[0::ny], color = c, linewidth = lthck)
    ax[1].plot(x / zscal, Y[1::ny], color = c, linewidth = lthck)
    ax[2].plot(x / zscal, Y[2::ny], color = c, linewidth = lthck)

    Y += dY * rmp

    Y[0::ny] = np.where(Y[0::ny] < 0, 0, Y[0::ny])
    Y[1::ny] = np.where(Y[1::ny] < T_floor/2, T_floor/2, Y[1::ny])

    plt.savefig('/tmp/fr/{:05d}.png'.format(it+1), dpi=160)
    plt.close(fig)
# from tempfile import mkstemp
# fid,fn = mkstemp('.svg')
#
# print 'saved as: {}'.format(fn)
# plt.show()
