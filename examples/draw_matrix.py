#!/usr/bin/env python
# coding: utf-8

from diskvert import *
import numpy as np

# initialize exapmle accretion disk
dv_init_disk(10,0.01,10)
dv_init_ss73(0.01)
dv_eval_globals()

# prepare memtory
nz = 2**10
z  = np.ndarray(nz)
y  = np.ndarray((4,nz),  order = 'F')
dy = np.ndarray((4,nz),  order = 'F')
a  = np.ndarray((8,nz),  order = 'F')

# run integration
dv_run_ss73(z,nz,y,dy,a)

# import ctypes
from ctypes import CDLL, POINTER, c_double, c_float, c_int
from numpy.ctypeslib import ndpointer

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
dv_init_ss73(0.03)
dv_eval_globals()

# draw
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from matplotlib.cm import get_cmap


fig,ax = plt.subplots(2, 2, figsize=(12,12))
cmap = get_cmap('Greens')

dv_generate_coefficients(x,nx,Y,ny,A,ny,M)
ax[0,0].set_title('Matrix')
ax[0,0].pcolor(M, norm = SymLogNorm(1, vmin = -1e2, vmax = 1e2), cmap='RdBu_r')
ax[0,0].set_xlabel('Y')
ax[0,0].set_xlim([0,len(Y)])
ax[0,0].set_ylabel('eq')
ax[0,0].set_ylim([0,len(A)])
ax[0,0].format_coord = lambda x,y : "x = %d, y = %d, z = %.3e" % (x, y, M[int(y),int(x)])

ax[0,1].set_title('rho')
ax[1,0].set_title('T')
ax[1,1].set_title('flux')

niter = 12
for it in range(niter):
    # generate the matrix and solve
    dv_generate_coefficients(x,nx,Y,ny,A,ny,M)
    dY = np.linalg.solve(M,-A)

    t = it / (niter - 1.0)
    ax[0,1].plot(x, Y[0::ny], color = cmap(t))
    ax[1,0].plot(x, Y[1::ny], color = cmap(t))
    ax[1,1].plot(x, Y[2::ny], color = cmap(t))

    Y = Y + dY * np.where( t < 0.75, t / 0.75, 1.0 )

    Y[1::ny] = np.where(Y[1::ny] < T_floor, T_floor, Y[1::ny])

from datetime import datetime
import tempfile
import os

fn = os.path.join(tempfile.mkdtemp(),
    '{}.png'.format(datetime.now().strftime('%Y%m%d%H%M%S')))
plt.savefig(fn)
print 'saved as: {}'.format(fn)
plt.show()
