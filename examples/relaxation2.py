#!/usr/bin/env python
# coding: utf-8

from diskvert import *
import numpy as np

# initialize exapmle accretion disk
dv_init_disk(10,0.001,10)
dv_init_ss73(0.1)
dv_eval_globals()

# prepare memtory
nz = 2**10
z  = np.ndarray(nz)
y  = np.ndarray((4,nz),  order = 'F')
dy = np.ndarray((4,nz),  order = 'F')
a  = np.ndarray((8,nz),  order = 'F')

# run integration
dv_run_ss73(z,nz,y,dy,a)

# interpolate initial values for new model
nx = 500
ny = 3
x = np.linspace(0, z.max(), nx)
Y = np.ndarray(nx*ny)
# Y[0::ny] = np.interp(x, z[::-1], a[2,::-1])
# Y[1::ny] = np.interp(x, z[::-1], a[1,::-1])
# Y[2::ny] = np.interp(x, z[::-1], y[2,::-1])
Y[0::ny] = np.linspace(0.02,1e-16,nx)
Y[1::ny] = np.linspace(3*a[1,0],a[1,0],nx)
Y[2::ny] = np.linspace(0,y[2,0],nx)

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


fig,ax = plt.subplots(1, 3, figsize=(21,7))
cmap = get_cmap('YlOrRd')

ax[0].set_title('rho')
ax[1].set_title('T')
ax[2].set_title('flux')


def ramp(it, niter, n = 1, ramp_0 = 0):
    t = np.clip(it / (niter - 1.0), 0, 1)
    return ramp_0 + (1 - ramp_0) * t**n

def ramp_smooth(it, niter, ramp0 = 0):
    t = np.clip(it / (niter - 1.0), 0, 1)
    return ramp0 + (1 - ramp0) * (3 * t**2 - 2 * t**3)

def ramp_trig(it, niter, ramp0 = 0):
    t = np.clip(it / (niter - 1.0), 0, 1)
    return ramp0 + (1 - ramp0) * (1 - np.cos(np.pi * t)) / 2

niter = 48

dY = np.ndarray(nx*ny)
special = [0, niter-1, 2*niter-1, 3*niter-1]
for it in range(niter * 3):
    # generate the matrix and solve
    dv_generate_corrections(x,nx,Y,dY,ny)

    rmp = [ ramp(it, niter, n, 0.01) for n in [2,3,1] ]

    dY = np.nan_to_num(dY)
    err = np.sqrt(np.sum( np.nan_to_num(dY/Y)**2 ))
    print "{:5d} {:8.2e} [{:s}]".format(it,err, ",".join([ "{:5.3f}".format(r) for r in rmp ]))

    err_check = err < 1e-12

    c = 'black' if it in special or err_check else cmap(ramp(it, niter))
    lthck = 1.25 if it in special or err_check else 0.8

    ax[0].plot(x, Y[0::ny], color = c, linewidth = lthck)
    ax[1].plot(x, Y[1::ny], color = c, linewidth = lthck)
    ax[2].plot(x, Y[2::ny], color = c, linewidth = lthck)

    if err_check:
        print 'precision reached!'
        break

    for i in range(ny):
        Y[i::ny] += dY[i::ny] * rmp[i]

    Y[0::ny] = np.where(Y[0::ny] < 0, 0, Y[0::ny])
    Y[1::ny] = np.where(Y[1::ny] < T_floor, T_floor, Y[1::ny])

from tempfile import mkstemp
fid,fn = mkstemp('.png')

plt.savefig(fn)
print 'saved as: {}'.format(fn)
plt.show()
