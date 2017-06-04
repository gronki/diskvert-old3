
# coding: utf-8

# In[43]:

from diskvert import *
import numpy as np

dv_init_disk(10,0.01,10)
dv_init_ss73(0.02)
dv_eval_globals()

nz = 1000
z = np.ndarray(nz, np.float64)
y = np.ndarray((nz,4), np.float64)
dy = np.ndarray((nz,4), np.float64)
a = np.ndarray((nz,8), np.float64)

dv_run_ss73(z,nz,y,dy,a)


# In[44]:

from ctypes import CDLL, POINTER, c_double, c_float, c_int
from numpy.ctypeslib import ndpointer


# In[45]:

dv_generate_coefficients.restype = None
dv_generate_coefficients.argtypes = [
    ndpointer(c_double,1),
    c_int,
    ndpointer(c_double,1),
    c_int,
    ndpointer(c_double,2),
    ndpointer(c_double,1),
]


# In[46]:

nx = 30
ny = 3
x = np.linspace(0, z.max(), nx)
Y = np.ndarray(nx*ny)
Y[0::ny] = np.interp(x[::-1],z[::-1],a[:,2])
Y[1::ny] = np.interp(x[::-1],z[::-1],a[:,1])
Y[2::ny] = np.interp(x[::-1],z[::-1],y[:,2])


# In[47]:

A = np.zeros((nx*ny))
M = np.zeros((nx*ny,nx*ny))


# In[48]:

#dv_init_ss73(0.04)
#dv_eval_globals()
dv_generate_coefficients(x,nx,Y,ny,M,A)


# In[49]:

import matplotlib.pyplot as plt
plt.pcolor(np.arcsinh(M * 1e8), vmin=-50, vmax=50, cmap='seismic')
plt.ylabel('Y')
plt.xlabel('eq')
plt.colorbar()
plt.show()


# In[ ]:
