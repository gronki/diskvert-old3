
# coding: utf-8

from diskvert import libdv
import numpy as np

libdv.init_disk(10,0.01,10)
libdv.init_ss73(0.1)
libdv.eval_globals()

nz = 1000
z = np.ndarray(nz, np.float64)
y = np.ndarray((nz,4), np.float64)
dy = np.ndarray((nz,4), np.float64)
a = np.ndarray((nz,8), np.float64)

libdv.run_ss73(z,nz,y,dy,a)

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
