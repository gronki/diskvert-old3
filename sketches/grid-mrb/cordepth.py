# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from diskvert import *
from sys import argv
from re import sub

for fn in argv[1:]:
    d,p = col2python(fn)
    z = d['z']
    T = d['temp']
    Tz = np.diff(T) / np.diff(z)
    zz = (z[1:] + z[:-1]) / 2

    zcor = 0
    for i in range(1,len(Tz)):
        if Tz[i-1] <= 0 and Tz[i] >= 0:
            zcor = (Tz[i]*zz[i-1] - Tz[i-1]*zz[i]) / (Tz[i] - Tz[i-1])
            break
    assert zcor != 0

    with open(sub(r'(\.dat|\.tar\.gz)$','.cor.txt',fn),'w') as f:
        strfmt = "{1:12s} {2:14.5e} # {0}\n"
        f.write(strfmt.format('height of the temperature inversion [cm]',
            'zcor', zcor))
        f.write(strfmt.format('height of the temperature inversion [H]',
            'hcor', zcor / p.zscale))
        taucor = np.interp(zcor,d['z'],d['tau'])
        f.write(strfmt.format('optical depth of the corona', 'taucor', taucor))
        fraddisk = np.interp(zcor,d['z'],d['frad'])
        f.write(strfmt.format('amout of radiation produced in the disk',
            'fraddisk', fraddisk))
        fradmax = np.max(d['frad'])
        chicor = 1 - fraddisk / fradmax
        f.write(strfmt.format('fraction of radiation produced in the corona',
            'chicor', chicor))
