from math import atan,sqrt

cgs_pi = 4*atan(1e0)
cgs_boltz = 1.3806581212e-16
cgs_c = 2.99792458e10
cgs_h = 6.62607554040e-27
cgs_graw = 6.67259858585e-8
cgs_mel = 9.10938975454e-28
cgs_hbar= 0.5 * cgs_h / cgs_pi
cgs_stef = 2 * cgs_pi**5 * cgs_boltz**4 / ( 15 * cgs_h**3 * cgs_c**2 )
cgs_a = 4 * cgs_stef / cgs_c
cgs_alpha = 1 / 137.036e0
cgs_qel = sqrt( cgs_alpha * cgs_hbar * cgs_c )
cgs_thomson = 8 * cgs_pi / 3 * cgs_qel**4 / ( cgs_c**4 * cgs_mel**2 )
cgs_elrad = cgs_qel**2 / ( cgs_c**2 * cgs_mel )
cgs_mhydr = 1.6733e-24
cgs_kapes1 = cgs_thomson / cgs_mhydr
cgs_k = cgs_boltz
cgs_k_over_mh = cgs_boltz / cgs_mhydr
cgs_mh_over_k = cgs_mhydr / cgs_boltz
cgs_k_over_mec = cgs_boltz / ( cgs_mel * cgs_c )
cgs_mec_over_k = ( cgs_mel * cgs_c ) / cgs_boltz
cgs_k_over_mec2 = cgs_boltz / ( cgs_mel * cgs_c**2 )
cgs_mec2_over_k = ( cgs_mel * cgs_c**2 ) / cgs_boltz
cgs_msun = 1.99e33
cgs_lsun = 3.9e33
cgs_kapes = 0.34e0
