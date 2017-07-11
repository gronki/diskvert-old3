#!/usr/bin/env python
# coding: utf-8

import numpy as np
from multiprocessing import cpu_count

ncpu = cpu_count()
radii = np.logspace(np.log10(3.16), np.log10(100), 12)

s0 = open('input.par','r').read()

with open('run.sh','w') as frun:
    frun.write("#!/usr/bin/env bash\n\n")
    for ir,r in zip(range(len(radii)), radii):
        ilab = "{0:03d}".format(ir+1)
        with open('{}.par'.format(ilab),'w') as f:
            f.write(s0 + "\nradius {:.4e}".format(r))
        frun.write("bash job.sh {0} &\n".format(ilab))
        if (ir+1) % ncpu == 0: frun.write("wait\n\n")
    frun.write("wait\n")
