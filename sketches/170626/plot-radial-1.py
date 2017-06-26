import matplotlib.pyplot as plt
import numpy as np

d = np.loadtxt('profile.dat')

fig, axes = plt.subplots(1, 3, figsize = (12,4))

ax = axes[0]
ax.plot(d[:,1], d[:,2], label = '3 zone')
ax.plot(d[:,1], d[:,3], label = 'iter')
ax.set_title('rho')

ax = axes[1]
ax.plot(d[:,1], d[:,4], label = '3 zone')
ax.plot(d[:,1], d[:,5], label = 'iter')
ax.set_title('T')

ax = axes[2]
ax.plot(d[:,1], d[:,6], label = '3 zone')
ax.plot(d[:,1], d[:,7], label = 'iter')
ax.set_title('H')

for ax in axes:
    ax.loglog()
    ax.legend(fontsize = 9)

plt.savefig('plot.png')
plt.show()
