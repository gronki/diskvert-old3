# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
from diskvert import *

dw,pw = col2python('MWFK.dat')
dz,pz = col2python('MZFK.dat')
dc,pc = col2python('MCFK.dat')

plt.yscale('log')
plt.plot(dc['h'],dc['trad'], label = 'trad', color = '#C0C0C0')
plt.plot(dc['h'],dc['temp'], label = 'corona', color = 'black', linewidth = 2)
plt.plot(dz['h'],dz['temp'], label='new compton', color = 'red')
plt.plot(dw['h'],dw['temp'], label='old compton', color = 'blue')
plt.legend()
plt.show()
