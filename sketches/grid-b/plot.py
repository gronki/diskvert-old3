from diskvert.pyminiconf import pyminiconf
from sys import argv
from matplotlib.pyplot import subplots, show, savefig
from numpy import ndarray

n = len(argv[1:])
betacor = ndarray(n)
tau_cor = ndarray(n)

for i,fn in zip(range(n),argv[1:]):
    p = pyminiconf(open(fn,'r'))
    betacor[i] = p.betacor
    tau_cor[i] = p.tau_cor
    print p.betacor

fig,ax = subplots()
ax.plot(betacor,tau_cor)
ax.set_xlabel('beta in temp minimum')
ax.set_ylabel('tau of the corona')
ax.loglog()
savefig('tau_vs_beta.png')
show()

