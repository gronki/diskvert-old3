
# coding: utf-8

from diskvert import *
import numpy as np

libdv.init_disk(10,0.01,10)
libdv.eval_globals()

nz = 100
z1 = np.ndarray(nz, np.float64)
z2 = np.ndarray(nz, np.float64)
z3 = np.ndarray(nz, np.float64)
hmax = 10

libdv.grid(GRID_LINEAR, hmax, z1, nz)
libdv.grid(GRID_LOG, hmax, z2, nz)
libdv.grid(GRID_ASINH, hmax, z3, nz)

import matplotlib.pyplot as plt

plt.figure(figsize = (20,10))

plt.subplot(1,2,1)
plt.title('grid locations')
plt.plot(z1, label='linear')
plt.plot(z2, label='log')
plt.plot(z3, label='asinh')
plt.legend(loc = 'upper left')

plt.subplot(1,2,2)
plt.title('grid delta')
plt.plot(z1[1:] - z1[0:-1], label='linear')
plt.plot(z2[1:] - z2[0:-1], label='log')
plt.plot(z3[1:] - z3[0:-1], label='asinh')
plt.legend(loc = 'upper left')

plt.show()
