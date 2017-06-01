from sympy import Symbol

alpha   = Symbol('alpha')
eta     = Symbol('eta')
Omega   = Symbol('Omega')
F_acc   = Symbol('F_acc')
mu      = Symbol('mu')
kappa_es = Symbol('kappa_es')
m_el     = Symbol('m_el')
m_H      = Symbol('m_H')
c        = Symbol('c')
sigma    = Symbol('sigma')
k        = Symbol('k_B')

global_variables = {
    Omega:  'Omega',
    F_acc:  'flux_acc',
    kappa_es:   'cgs_kapes',
    m_el:       'cgs_mel',
    m_H:        'cgs_mhydr',
    c:          'cgs_c',
    sigma:      'cgs_stef',
    k:          'cgs_boltz',
    mu:         'miu',
}
