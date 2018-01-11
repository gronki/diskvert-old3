import numpy as np

nrad = 24
N2 = 24

mbh = 10
mdot = 0.01
mdots = np.logspace(-4, -1, N2)
radii = np.logspace(0.5, 2.4, nrad)

alpha = 0.01

eta = 0.5
etas = np.logspace(np.log10(alpha), 0, N2)
betas = 2 * etas / alpha - 1

fn = lambda i,j: "R{:02d}X{:02d}".format(i+1, j+1)

def DPA(D, k, mask = True):
    from numpy import vectorize
    d0 = vectorize(lambda d: getattr(d,k) if d.converged else np.nan)(D)
    if not mask: return d0
    from numpy.ma import masked_where
    dm = vectorize(lambda d: not d.converged)(D)
    return masked_where(dm, d0)

def read2Ddataset(nx, ny, fnpath):
    from diskvert import pyminiconf
    from numpy import ndarray
    D = ndarray((nx,ny), dtype = object)
    for i in range(nx):
        for j in range(ny):
            with open(fnpath(i,j),'r') as f:
                D[i,j] = pyminiconf(f)
    return D
    
if __name__ == '__main__':
    par = """
    mbh     {mbh}
    mdot    {mdot}
    alpha   {alpha}
    radius  {radius}
    zeta    {eta}
    """

    fnlist = open('jobs.lst','w')

    for i in range(nrad):
        for j in range(N2):
            with open('par/' + fn(i,j) + '.par','w') as f:
                f.write(par.format(
                    mbh = mbh,
                    mdot = mdot,
                    alpha = alpha,
                    radius = radii[i],
                    eta = etas[j]
                ))
                fnlist.write(fn(i,j) + '\n')

    fnlist.close()
