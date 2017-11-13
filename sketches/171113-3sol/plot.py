from diskvert import col2python
import matplotlib.pyplot as plt
from numpy import exp

# relaksacja
dr,pr = col2python('relx.dat')
# runge kutta + multi temp
dm,pm = col2python('rkm.dat')

plt.figure(figsize = (12,9))

plt.plot(dr['h'], dr['trad'], color = '#CBCBCB', linewidth = 4, label = 'Trad (relax)')
plt.plot(dr['h'], dr['temp'], color = '#8B8B8B', linewidth = 4, label = 'Tgas (relax)')

plt.plot(dm['h'], exp(dm['tbil1']), color = '#F7DE36', linewidth = 1.5, label = 'Tgas1 (RK4)')
plt.plot(dm['h'], exp(dm['tbil2']), color = '#F76B26', linewidth = 1.5, label = 'Tgas2 (RK4)')
plt.plot(dm['h'], exp(dm['tbil3']), color = '#E40059', linewidth = 1.5, label = 'Tgas3 (RK4)')

plt.xlabel('z / H')
plt.ylabel('T [K]')

plt.yscale('log')
plt.ylim(1e6,1e9)

plt.legend()
plt.savefig('temp.png')
