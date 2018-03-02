import numpy as np

sym_labels = ['ZNVF', '$\\beta 5$', '$\\beta 4$', '$\\beta 3$', '$\\beta 2$', '$\\beta 1$']

sym_betas = np.array([70, 40, 17, 2, 0.7, 0.6])
sym_betas_err = np.array([20, 10, 14, 1, 0.7, 0.6])

sym_alphas = np.array([0.016, 0.029, 0.09, 0.30, 0.30, 0.4])
sym_alphas_err = np.array([0.002, 0.003, 0.01, 0.10, 0.20, 0.30])

sym_etas = np.array([ 0.029, 0.043, 0.068, 0.12, 0.21, 0.30 ])
sym_zetas = np.array([ 0.5, 0.6, 0.7, 0.4, 0.08, 0.00 ])

sym_nus = 2 * sym_zetas / sym_alphas
# sym_nus_err = sym_nus * (sym_alphas_err / sym_alphas + sym_zetas_err / sym_zetas)
