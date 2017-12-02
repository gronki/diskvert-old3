
from numpy import linspace, logspace, log10

mbh = 10
mdot = 0.01
alpha = 0.007
radius = 9.22

izip = lambda x: zip(range(len(x)), x)

for inu, nu in izip(logspace(log10(0.1), log10(20), 12)):
    for izet, zeta in izip(logspace(log10(alpha), 0.0, 12)):
        beta = 2 * zeta / alpha - 1 + nu
        if beta <= 0: continue
        fn = 'B{:02d}N{:02d}'.format(izet + 1, inu + 1)
        with open('par/' + fn, 'w') as f:
            f.write(    "# model    {}\n" \
                        "  mbh      {}\n" \
                        "  mdot     {}\n" \
                        "  radius   {}\n" \
                        "  alpha    {}\n" \
                        "  zeta     {} # nr {}\n" \
                        "# beta     {}\n" \
                        "  nu       {} # nr {}\n" \
            .format(fn, mbh, mdot, radius, alpha, zeta, izet + 1, beta, nu, inu + 1))
