from factory import *
import os.path

__libdv = CDLL('libdiskvert.so')

GRID_LINEAR = 1 << 1
GRID_LOG    = 2 << 1
GRID_ASINH  = 3 << 1
GRID_REVERSE = 1

EQUATION_EQUILIBR   = 0
EQUATION_COMPTON    = 1
EQUATION_BALANCE    = 2

dv_init_disk = CFACTORY(__libdv,'dv_init_disk', [
    (c_double, 'in', 'mass'),
    (c_double, 'in', 'accretion'),
    (c_double, 'in', 'radius'),
])

dv_init_abun = CFACTORY(__libdv,'dv_init_abun', [
    (c_double, 'in', 'X'),
    (c_double, 'in', 'Z'),
])

dv_eval_globals = CFACTORY(__libdv,'dv_eval_globals')

# dv_grid = CFACTORY(__libdv,'grid', [
#     (c_int, 'in', 'type'),
#     (c_double, 'in', 'range'),
#     (ndpointer(c_double,1), 'inout', 'values'),
#     (c_int, 'in', 'N')
# ])

dv_init_ss73 = CFACTORY(__libdv,'init_ss73', [
    (c_double, 'in', 'alpha')
])

dv_run_ss73 = CFACTORY(__libdv,'run_ss73', [
    (ndpointer(c_double,1), 'inout', 'z'),
    (c_int, 'in', 'nz'),
    (ndpointer(c_double,2), 'inout', 'Y'),
    (ndpointer(c_double,2), 'inout', 'dY'),
    (ndpointer(c_double,2), 'inout', 'A'),
])

dv_generate_coefficients = CFACTORY(__libdv,'generate_coefficients', [
    (ndpointer(c_double,1), 'in', 'z'),
    (c_int, 'in', 'nz'),
    (ndpointer(c_double,1), 'in', 'Y'),
    (c_int, 'in', 'ny'),
    (ndpointer(c_double,1), 'inout', 'A'),
    (c_int, 'in', 'na'),
    (ndpointer(c_double,2), 'inout', 'M'),
])

dv_generate_corrections = CFACTORY(__libdv,'generate_corrections', [
    (ndpointer(c_double,1), 'in', 'z'),
    (c_int, 'in', 'nz'),
    (ndpointer(c_double,1), 'in', 'Y'),
    (ndpointer(c_double,1), 'inout', 'dY'),
    (c_int, 'in', 'ny'),
])

dv_query_f = CFACTORY(__libdv,'dv_query_f', [
    (POINTER(c_char), 'in', 'key'),
], c_float)
