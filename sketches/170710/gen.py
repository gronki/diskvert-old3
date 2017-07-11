#!/usr/bin/env python
# coding: utf-8

import numpy as np
from multiprocessing import cpu_count
from diskvert.pyminiconf import pyminiconf

ncpu = 12 # cpu_count()

p = pyminiconf(open('input.par','r'))

beta_max = 2.0 / p.alpha - 1.0
beta = np.logspace(np.log10(beta_max), np.log10(0.2), 24)
X = p.alpha * (beta + 1) / 2

s0 = open('input.par','r').read()

with open('run.sh','w') as frun:
    frun.write("#!/usr/bin/env bash\n\n")
    for ix,x in zip(range(len(X)), X):
        ilab = "{0:04d}".format(ix+1)
        with open('D.{}.par'.format(ilab),'w') as f:
            f.write(s0 + "\nzeta {:.4e}".format(x))
        frun.write("bash job.sh {0} &\n".format(ilab))
        if (ix+1) % ncpu == 0: frun.write("wait\n\n")
    frun.write("wait\n")
