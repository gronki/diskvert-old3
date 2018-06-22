import numpy as np

def random_magpars():
    eta = 10**np.random.uniform(-2, -0.3)
    beta = 10**np.random.uniform(0, 2)
    alpha = 2 * eta / (beta + 1)
    nu = 10**np.random.uniform(-1,1)
    return alpha, eta, nu

if __name__ == '__main__':
    alpha, eta, nu = random_magpars()
    print u'alpha, eta, nu = {:.5f}, {:.3f}, {:.3f}'.format(alpha, eta, nu)
