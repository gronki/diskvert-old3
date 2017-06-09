#!/usr/bin/env python
# coding: utf-8

from diskvert import *
import numpy as np

# initialize exapmle accretion disk
dv_init_disk(10,0.01,10)
dv_init_ss73(0.2)
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
nx = 2**8
ny = 3
x = np.linspace(0, z.max(), nx)
Y = np.ndarray(nx*ny)
Y[0::ny] = np.interp(x, z[::-1], a[2,::-1])
Y[1::ny] = np.interp(x, z[::-1], a[1,::-1])
Y[2::ny] = np.interp(x, z[::-1], y[2,::-1])

T_floor = a[1,0]

# matrix
A = np.zeros((nx*ny))
M = np.zeros((nx*ny,nx*ny), order = 'F')

# change parameters a little bit
dv_init_ss73(0.02)
dv_eval_globals()

# draw
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from matplotlib.cm import get_cmap


fig,ax = plt.subplots(1, 3, figsize=(21,7))
cmap = get_cmap('Oranges')

ax[0].set_title('rho')
ax[1].set_title('T')
ax[2].set_title('flux')

def ramp1(it, niter, ramp_0 = 0):
    t = np.clip(it / (niter - 1.0), 0, 1)
    return ramp_0 + (1 - ramp_0) * t
def ramp2(it, niter, ramp_0 = 0):
    t = np.clip(it / (niter - 1.0), 0, 1)
    return ramp_0 + (1 - ramp_0) * t**2
def ramp3(it, niter, ramp_0 = 0):
    t = np.clip(it / (niter - 1.0), 0, 1)
    return ramp_0 + (1 - ramp_0) * t**3
def ramp_poly(it, niter, ramp0 = 0):
    t = np.clip(it / (niter - 1.0), 0, 1)
    return ramp0 + (1 - ramp0) * (3 * t**2 - 2 * t**3)
def ramp_trig(it, niter, ramp0 = 0):
    t = np.clip(it / (niter - 1.0), 0, 1)
    return ramp0 + (1 - ramp0) * (1 - np.cos(np.pi * t)) / 2

niter = 12
niter_full = 3

for it in range(niter + niter_full - 1):
    # generate the matrix and solve
    dv_generate_coefficients(x,nx,Y,ny,A,ny,M)
    dY = np.linalg.solve(M,-A)

    rmp = ramp2(it, niter, 0.01)

    c = 'black' if it == 0 or it == niter-1 else cmap(rmp)

    ax[0].plot(x, Y[0::ny], color = c)
    ax[1].plot(x, Y[1::ny], color = c)
    ax[2].plot(x, Y[2::ny], color = c)

    err = np.sqrt(np.sum( np.nan_to_num(dY/Y)**2 ))
    print "{:5d} {:8.2e} {:5.2f} {:8.2e}".format(it,err,rmp, err*rmp)

    Y = Y + dY * rmp

    Y[0::ny] = np.where(Y[0::ny] < 0, 0, Y[0::ny])
    Y[1::ny] = np.where(Y[1::ny] < T_floor, T_floor, Y[1::ny])

from tempfile import mkstemp
fid,fn = mkstemp('.png')

plt.savefig(fn)
print 'saved as: {}'.format(fn)
plt.show()
