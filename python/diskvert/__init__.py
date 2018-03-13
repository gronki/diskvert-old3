# coding: utf-8

from diskvert.cgs import *
from diskvert.col2python import col2python
from diskvert.pyminiconf import pyminiconf, pyminiconf_dict
import diskvert.cm
from plotrx import rxread

zeta = lambda alfa,beta: alfa*(1 + beta)/2
