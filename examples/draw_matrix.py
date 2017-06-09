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

# interpolate initial values for new model
nx = 8
ny = 3
x = np.linspace(0, z.max(), nx)
Y = np.ndarray(nx*ny)
Y[0::ny] = np.interp(x, z[::-1], a[2,::-1])
Y[1::ny] = np.interp(x, z[::-1], a[1,::-1])
Y[2::ny] = np.interp(x, z[::-1], y[2,::-1])

# matrix
A = np.zeros((nx*ny))
M = np.zeros((nx*ny,nx*ny), order = 'F')

# draw
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from matplotlib.cm import get_cmap


fig,ax = plt.subplots(figsize=(12,12))
cmap = get_cmap('Greens')

dv_generate_coefficients(x,nx,Y,ny,A,ny,M)
ax.set_title('Matrix')
ax.pcolor(M, norm = SymLogNorm(1, vmin = -1e2, vmax = 1e2), cmap='RdBu_r')
ax.set_xlabel('Y')
ax.set_xlim([0,len(Y)])
ax.set_ylabel('eq')
ax.set_ylim([0,len(A)])
ax.format_coord = lambda x,y : "x = %d, y = %d, z = %.3e" % (x, y, M[int(y),int(x)])

from tempfile import mkstemp
fid,fn = mkstemp('.png')

plt.savefig(fn)
print 'saved as: {}'.format(fn)
plt.show()
