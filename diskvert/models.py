from ctypes import CDLL, CFUNCTYPE, POINTER, c_double, c_float, c_int
from numpy.ctypeslib import ndpointer
import os.path

__libdv = CDLL('libdiskvert.so')

def FACTORY(name, args = None, restype = None):
    flags = { 'in':1, 'out':2, 'inout':3 }
    args = args if args else []
    argtypes = tuple( a_type for a_type,a_intent,a_name in args )
    arglist = tuple( (flags[a_intent], a_name) for a_type,a_intent,a_name in args )
    proto = CFUNCTYPE(restype, *argtypes)
    f = proto((name,__libdv), arglist)
    f.__doc__ = '{}({})'.format(name, ', '.join([
        '{} {}'.format(a_type.__name__,a_name) \
        for a_type,a_intent,a_name in args  ]))
    return f

GRID_LINEAR = 1 << 1
GRID_LOG    = 2 << 1
GRID_ASINH  = 3 << 1
GRID_REVERSE = 1

EQUATION_EQUILIBR   = 0
EQUATION_COMPTON    = 1
EQUATION_BALANCE    = 2

dv_init_disk = FACTORY('init_disk', [
    (c_double, 'in', 'mass'),
    (c_double, 'in', 'accretion'),
    (c_double, 'in', 'radius'),
])

dv_init_abun = FACTORY('init_abun', [
    (c_double, 'in', 'X'),
    (c_double, 'in', 'Z'),
])

dv_eval_globals = FACTORY('eval_globals')

dv_grid = FACTORY('grid', [
    (c_int, 'in', 'type'),
    (c_double, 'in', 'range'),
    (ndpointer(c_double,1), 'inout', 'values'),
    (c_int, 'in', 'N')
])

dv_init_ss73 = FACTORY('init_ss73', [
    (c_double, 'in', 'alpha')
])

dv_run_ss73 = FACTORY('run_ss73', [
    (ndpointer(c_double,1), 'inout', 'z'),
    (c_int, 'in', 'nz'),
    (ndpointer(c_double,2), 'inout', 'y'),
    (ndpointer(c_double,2), 'inout', 'dy'),
    (ndpointer(c_double,2), 'inout', 'a'),
])

dv_generate_coefficients = FACTORY('generate_coefficients', [
    (ndpointer(c_double,1), 'in', 'z'),
    (c_int, 'in', 'nz'),
    (ndpointer(c_double,1), 'in', 'y'),
    (c_int, 'in', 'ny'),
    (ndpointer(c_double,2), 'inout', 'M'),
    (ndpointer(c_double,1), 'inout', 'A'),
])
