from ctypes import CDLL, POINTER, c_double, c_float, c_int
from numpy.ctypeslib import ndpointer
import os.path

libdv = CDLL('libdiskvert.so')

GRID_LINEAR = 1 << 1
GRID_LOG    = 2 << 1
GRID_ASINH  = 3 << 1
GRID_REVERSE = 1
GRID_REVERSE_MASK = 1
GRID_TYPE_MASK = ~GRID_REVERSE_MASK

EQUATION_EQUILIBR   = 0
EQUATION_COMPTON    = 1
EQUATION_BALANCE    = 2

libdv.init_disk.restype = None
libdv.init_disk.argtypes = [
    c_double,
    c_double,
    c_double,
]

libdv.init_abun.restype = None
libdv.init_abun.argtypes = [
    c_double,
    c_double,
]

libdv.eval_globals.restype = None
libdv.eval_globals.argtypes = None

libdv.grid.restype = None
libdv.grid.argtypes = [
    c_int,
    c_double,
    ndpointer(c_double,1),
    c_int,
]

libdv.init_ss73.restype = None
libdv.init_ss73.argtypes = [ c_double ]

libdv.run_ss73.restype = None
libdv.run_ss73.argtypes = [
    ndpointer(c_double,1),
    c_int,
    ndpointer(c_double, 2),
    ndpointer(c_double, 2),
    ndpointer(c_double, 2),
]

# libdv.single_ss73
#
# libdv.init_m1
# libdv.run_m1
# libdv.single_m1
