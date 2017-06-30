from factory import *

__libdv = CDLL('libdiskvert.so')

alf_init = CFACTORY(__libdv,'init_alpha', [
    (POINTER(c_double), 'in', 'mbh'),
    (POINTER(c_double), 'in', 'mdot'),
    (POINTER(c_double), 'in', 'r'),
    (POINTER(c_double), 'in', 'alpha')
])

alf_run = CFACTORY(__libdv,'run_alpha', [
    (ndpointer(c_double,1), 'inout', 'z'),
    (POINTER(c_int), 'in', 'nz'),
    (ndpointer(c_double,2), 'inout', 'Y'),
    (ndpointer(c_double,2), 'inout', 'dY'),
    (ndpointer(c_double,2), 'inout', 'A'),
])

tempmethod = CFACTORY(__libdv,'tempmethod', [
    (c_char, 'in', 'method'),
])

alf_getn = CFACTORY(__libdv, 'alf_getn', [
    (POINTER(c_int), 'out', 'ny'),
    (POINTER(c_int), 'out', 'na'),
])

apxdisk2d = CFACTORY(__libdv, 'apxdisk2d', [
    (c_double, 'in', 'mass_black_hole'),
    (c_double, 'in', 'accretion_rate_mdot'),
    (ndpointer(c_double), 'in', 'radii'),
    (c_double, 'in', 'alpha'),
    (ndpointer(c_double), 'in', 'heights'),
    (ndpointer(c_double,2), 'inout', 'rho'),
    (ndpointer(c_double,2), 'inout', 'T'),
    (c_int, 'in', 'n_radii'),
    (c_int, 'in', 'n_heights'),
])

mrx_model = CFACTORY(__libdv, 'mrx_model', [
    (c_int, 'in', 'compton'),
    (c_int, 'in', 'magnetic'),
    (c_int, 'in', 'conduction'),
    (POINTER(c_int), 'out', 'nr'),
])

mrx_init = CFACTORY(__libdv, 'mrx_init', [
    (c_double, 'in', 'mbh'),
    (c_double, 'in', 'mdot'),
    (c_double, 'in', 'r'),
    (c_double, 'in', 'alpha'),
    (c_double, 'in', 'zeta'),
])

mrx_get_ny = CFACTORY(__libdv, 'mrx_get_ny', [
    (c_int, 'in', 'nr'),
    (POINTER(c_int), 'out', 'ny'),
])

mrx_matrix = CFACTORY(__libdv, 'mrx_matrix', [
    (c_int, 'in', 'nr'),
    (ndpointer(c_double,1), 'in', 'z'),
    (c_int, 'in', 'nz'),
    (ndpointer(c_double,1), 'in', 'Y'),
    (c_int, 'in', 'ny'),
    (ndpointer(c_double,1), 'inout', 'A'),
    (c_int, 'in', 'na'),
    (ndpointer(c_double,2), 'inout', 'M'),
])

mrx_advance = CFACTORY(__libdv, 'mrx_advance', [
    (c_int, 'in', 'nr'),
    (ndpointer(c_double,1), 'in', 'z'),
    (c_int, 'in', 'nz'),
    (ndpointer(c_double,1), 'in', 'Y'),
    (ndpointer(c_double,1), 'inout', 'dY'),
    (c_int, 'in', 'ny'),
])
