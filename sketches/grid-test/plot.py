#!/usr/bin/env python
# coding=utf-8

#------------------------------------------------------------------------------#

from diskvert import col2python
import numpy as np

#------------------------------------------------------------------------------#

lin060d, lin060p = col2python('data/lin060.dat')
lin120d, lin120p = col2python('data/lin120.dat')
lin240d, lin240p = col2python('data/lin240.dat')
lin480d, lin480p = col2python('data/lin480.dat')

lind, linp = col2python('data/lin.dat')
logd, logp = col2python('data/log.dat')
ashd, ashp = col2python('data/ash.dat')
gsqd, gsqp = col2python('data/gsq.dat')

log060d, log060p = col2python('data/log060.dat')
log120d, log120p = col2python('data/log120.dat')
log240d, log240p = col2python('data/log240.dat')
log480d, log480p = col2python('data/log480.dat')

#------------------------------------------------------------------------------#

import matplotlib.pyplot as plt

plt.figure(figsize = (10,10))

plt.subplot(221)
plt.plot(lin480d['h'], lin480d['frad'], label = u'rad (H = 480)')
plt.plot(lin480d['h'], lin480d['fmag'], '--', label = u'mag (H = 480)')
plt.plot(lin240d['h'], lin240d['frad'], label = u'rad (H = 240)')
plt.plot(lin240d['h'], lin240d['fmag'], '--', label = u'mag (H = 240)')
plt.plot(lin120d['h'], lin120d['frad'], label = u'rad (H = 120)')
plt.plot(lin120d['h'], lin120d['fmag'], '--', label = u'mag (H = 120)')
plt.plot(lin060d['h'], lin060d['frad'], label = u'rad (H = 060)')
plt.plot(lin060d['h'], lin060d['fmag'], '--', label = u'mag (H = 060)')
plt.xscale('log')
plt.xlim(0.1,None)
plt.legend(loc = 'best')

plt.subplot(223)
plt.plot(lin480d['h'], lin480d['temp'], label = u'linear (n = 480)')
plt.plot(lin240d['h'], lin240d['temp'], label = u'linear (n = 240)')
plt.plot(lin120d['h'], lin120d['temp'], label = u'linear (n = 120)')
plt.plot(lin060d['h'], lin060d['temp'], label = u'linear (n = 060)')
plt.yscale('log')
plt.xscale('log')
plt.xlim(0.1,None)
plt.legend(loc = 'best')

plt.subplot(222)
plt.plot(lind['i'], lind['temp'], label = u'linear')
plt.plot(logd['i'], logd['temp'], label = u'log')
plt.plot(ashd['i'], ashd['temp'], label = u'asinh')
plt.plot(gsqd['i'], gsqd['temp'], label = u'squared')
plt.yscale('log')
plt.legend(loc = 'best')

plt.subplot(224)
plt.plot(lin480d['h'], lin480d['rho'], label = u'linear (n = 480)')
plt.plot(lin240d['h'], lin240d['rho'], label = u'linear (n = 240)')
plt.plot(lin120d['h'], lin120d['rho'], label = u'linear (n = 120)')
plt.plot(lin060d['h'], lin060d['rho'], label = u'linear (n = 060)')
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.1,None)

plt.legend(loc = 'best')

#------------------------------------------------------------------------------#

clr = [
    '#F14629',
    '#48D53E',
    '#1B6BBF',
    '#B60AA0',
]

plt.figure(figsize = (10,10))

plt.subplot(221)
plt.plot(log480d['h'], log480d['frad'], color = clr[3], label = u'rad (H = 480)')
plt.plot(log480d['h'], log480d['fmag'], '--', color = clr[3], label = u'mag (H = 480)')
plt.plot(log240d['h'], log240d['frad'], color = clr[2], label = u'rad (H = 240)')
plt.plot(log240d['h'], log240d['fmag'], '--', color = clr[2], label = u'mag (H = 240)')
plt.plot(log120d['h'], log120d['frad'], color = clr[1], label = u'rad (H = 120)')
plt.plot(log120d['h'], log120d['fmag'], '--', color = clr[1], label = u'mag (H = 120)')
plt.plot(log060d['h'], log060d['frad'], color = clr[0], label = u'rad (H = 060)')
plt.plot(log060d['h'], log060d['fmag'], '--', color = clr[0], label = u'mag (H = 060)')
plt.xscale('log')
plt.xlim(0.1,None)
plt.legend(loc = 'best')

plt.subplot(223)
plt.plot(log480d['h'], log480d['temp'], label = u'linear (n = 480)')
plt.plot(log240d['h'], log240d['temp'], label = u'linear (n = 240)')
plt.plot(log120d['h'], log120d['temp'], label = u'linear (n = 120)')
plt.plot(log060d['h'], log060d['temp'], label = u'linear (n = 060)')
plt.yscale('log')
plt.xscale('log')
plt.xlim(0.1,None)
plt.legend(loc = 'best')

plt.subplot(222)
plt.plot(log480d['i'], log480d['temp'], label = u'480')
plt.plot(log240d['i'], log240d['temp'], label = u'240')
plt.plot(log120d['i'], log120d['temp'], label = u'120')
plt.plot(log060d['i'], log060d['temp'], label = u'60')
plt.yscale('log')
plt.legend(loc = 'best')

plt.subplot(224)
plt.plot(log480d['i'], log480d['rho'], label = u'480')
plt.plot(log240d['i'], log240d['rho'], label = u'240')
plt.plot(log120d['i'], log120d['rho'], label = u'120')
plt.plot(log060d['i'], log060d['rho'], label = u'60')
plt.yscale('log')
plt.legend(loc = 'best')

#------------------------------------------------------------------------------#

plt.show()

#------------------------------------------------------------------------------#
