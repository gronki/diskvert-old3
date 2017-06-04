
# coding: utf-8

from diskvert import *
import numpy as np

dv_init_disk(10,0.01,10)
dv_init_ss73(0.1)
dv_eval_globals()

nz = 1000
z = np.ndarray(nz, np.float64)
y = np.ndarray((nz,4), np.float64)
dy = np.ndarray((nz,4), np.float64)
a = np.ndarray((nz,8), np.float64)

dv_run_ss73(z,nz,y,dy,a)

import matplotlib.pyplot as plt

plt.subplot(121)
plt.plot(z,a[:,3],label='$T_{rad}$')
plt.plot(z,a[:,1],label='$T_{gas}$')
plt.yscale('log')
plt.legend()

plt.subplot(122)
plt.plot(z,a[:,2],label='$\\rho$')
plt.yscale('log')
plt.legend()

plt.show()
