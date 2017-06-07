#!/usr/bin/env python
# coding: utf-8

from diskvert import *
import numpy as np

# initialize exapmle accretion disk
dv_init_disk(10,0.01,10)
dv_init_ss73(0.02)
dv_eval_globals()

# prepare memtory
nz = 1000
z = np.ndarray(nz, np.float64)
y = np.ndarray((nz,4), np.float64)
dy = np.ndarray((nz,4), np.float64)
a = np.ndarray((nz,8), np.float64)

# run integration
dv_run_ss73(z,nz,y,dy,a)

# import ctypes
from ctypes import CDLL, POINTER, c_double, c_float, c_int
from numpy.ctypeslib import ndpointer

# interpolate initial values for new model
nx = 30
ny = 3
x = np.linspace(0, z.max(), nx)
Y = np.ndarray(nx*ny)
Y[0::ny] = np.interp(x[::-1],z[::-1],a[:,2])
Y[1::ny] = np.interp(x[::-1],z[::-1],a[:,1])
Y[2::ny] = np.interp(x[::-1],z[::-1],y[:,2])

# matrix
A = np.zeros((nx*ny))
M = np.zeros((nx*ny,nx*ny))

# change parameters a little bit
#dv_init_ss73(0.04)
#dv_eval_globals()

# generate the matrix
dv_generate_coefficients(x,nx,Y,ny,A,ny,M)

# draw
import matplotlib.pyplot as plt
plt.pcolor(np.arcsinh(M * 1e8), vmin=-50, vmax=50, cmap='seismic')
plt.ylabel('Y')
plt.xlabel('eq')
plt.colorbar()
plt.show()
