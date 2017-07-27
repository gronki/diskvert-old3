# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt

tx = np.linspace(0, 1, 1000)
# t = tx / (1 - tx)
t = (1 + tx) / (1 - tx)
t0 = (1 + t) * (1 + t**2)

fig = plt.figure(figsize = (12,9), dpi = 118)

plt.plot(tx, 4 / t0, label = '$4$', color = '#F1230E')
plt.plot(tx, (6*t-2) / t0, label = '$6t-2$', color = '#FFEC45')
plt.plot(tx, (1+3*t**2) / t0, label = '$1+3t^2$', color = '#E5FF45')
plt.plot(tx, t**3 / t0, label = '$t^3$', color = '#0CCD1A')

plt.plot(tx, (3+t**3) / t0, label = '$3+t^3$', color = '#49F6E2')
plt.plot(tx, t*(3+t**2) / t0, label = '$t(3+t^2)$', color = '#167DDD')
plt.plot(tx, (3*(1+t**2)/2+t**3) / t0, label = '$\\frac{3}{2}(1+t^2)+t^3$', color = '#500BC1')

plt.plot(tx, (1+t)*(1+t**2) / t0, color = 'black', linewidth = 1.6, label = '$1 + t + t^2 + t^3$')
plt.legend(loc = 'best')
plt.savefig('krzywe.png')
plt.show()
