
# coding: utf-8

# In[14]:

from diskvert import *
d,p = col2python('MCFX.dat')


# In[15]:

clr_blue = '#2342A7'
clr_green = '#07A82B'
clr_orange = '#F78517'
clr_orange_2 = '#834000'


# In[16]:

from numpy import ndarray
def integratetau(z,rho,kap,y):
    n = y.shape[0]
    no = ndarray(n)
    de = ndarray(n)
    de[n-1] = rho[n-1] * kap[n-1] * (z[n-1] - z[n-2])
    no[n-1] = de[n-1] * y[n-1]
    for i in range(1,n):
        dz = z[n - i] - z[n - (i+1)]
        rho_m = (rho[n - i] + rho[n - (i+1)]) / 2
        kap_m = (kap[n - i] + kap[n - (i+1)]) / 2
        y_m = (y[n - i] + y[n - (i+1)]) / 2
        no[n - (i+1)] = no[n - i] + dz * rho_m * kap_m * y_m
        de[n - (i+1)] = de[n - i] + dz * rho_m * kap_m
    return no / de


# In[20]:

from numpy import sqrt
from matplotlib.pyplot import subplots, show, savefig
fig,ax = subplots(figsize = (12,9))
# plot temperatures
ax.plot(d['h'], d['trad'], color = '#AAAAAA', linewidth = 1, label = u'Trad')
ax.plot(d['h'], d['temp'], color = 'black', linewidth = 2, label = u'Tgas')
# thermal kappa
ktot = d['ksct'] + d['kabs']
kthm = sqrt(d['kabs'] * ktot)
# plot average temperatures
ax.plot(d['h'], integratetau(d['z'], d['rho'], ktot, d['temp']), color = clr_orange, label = u'Tavg (over tau_tot)')
ax.plot(d['h'], integratetau(d['z'], d['rho'], kthm, d['temp']), color = clr_green, label = u'Tavg (over tau_th)')
ax.plot(d['h'], d['tavg'], linestyle = '--', color = clr_orange_2, label = u'Tavg (over tau, in Fortran)')
ax.axvline(p.hcor, linestyle = ':', color = '#D70A0A', label = 'temperature minimum')
ax.axvline(p.zphot / p.zscale, linestyle = '--', color = '#A2B1F6', label = 'photosphere (tau = 1)')
ax.axvline(p.ztherm / p.zscale, linestyle = '--', color = '#13A71E', label = 'thermalization depth')
# set limits
ax.set_xlim(0,30)
ax.set_ylim(1e6,1e8)
ax.set_yscale('log')
ax.set_xlabel('z / H')
ax.set_ylabel('temperature')
ax.legend(loc = 'best')
# save and show
savefig('tempav.png')
show()
