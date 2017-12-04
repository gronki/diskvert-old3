from diskvert.pyminiconf import pyminiconf
from sys import argv
from matplotlib.pyplot import subplots, show, savefig
from numpy import array, ndarray, nan

#-------------------------------------------------------------------------------

n = len(argv[1:])
l = []
for i,fn in zip(range(n),argv[1:]):
    l.append(pyminiconf(open(fn,'r')))

#-------------------------------------------------------------------------------

def f(o,a):
    from numpy import nan
    if hasattr(o,a):
        if type(getattr(o,a)) != str:
            return getattr(o,a)
    return nan

betacor = array([f(p,'betacor') for p in l])
tau_cor = array([f(p,'tau_cor') for p in l])
tcor = array([f(p,'tcor') for p in l])
chicor = array([f(p,'chicor') for p in l])

#-------------------------------------------------------------------------------

fig,ax = subplots()
ax.plot(betacor,tau_cor)

ax.set_xlabel('beta in temp minimum')
ax.set_xscale('log')
ax.set_xlim(1e-4,10)

ax.set_ylabel('optical depth')
ax.set_ylim(0,20)
savefig('tau_vs_beta.png')

#-------------------------------------------------------------------------------

fig,ax = subplots()
ax.plot(betacor,chicor)

ax.set_xlabel('beta in temp minimum')
ax.set_xscale('log')
ax.set_xlim(1e-4,10)

ax.set_ylabel('energy in corona')
ax.set_ylim(0,1)

savefig('chicor_vs_beta.png')

#-------------------------------------------------------------------------------

fig,ax = subplots()
ax.plot(betacor,chicor)

ax.set_xlabel('beta in temp minimum')
ax.set_xscale('log')
ax.set_xlim(1e-4,10)

ax.set_ylabel('corona temperature [keV]')
ax.set_ylim(0,2)

savefig('tcor_vs_beta.png')

#-------------------------------------------------------------------------------

fig,ax = subplots()
ax.plot(tau_cor, tcor / 1.16e7)

ax.set_xlabel('optical depth')
ax.set_xlim(0,20)

ax.set_ylabel('corona temperature [keV]')
ax.set_ylim(0,2)

savefig('tau_vs_tcor.png')

#-------------------------------------------------------------------------------

# show()
