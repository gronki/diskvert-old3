rho1 = 2.3049d-6*1.0/mbh*r**1.5d0/(alpha*cgs_kapes**3*mdot**2*f(r)**2)
T1 = 4.6189d+7*alpha**(-0.25d0)*cgs_kapes**(-0.25d0)*mbh**(-0.25d0)*r**( &
      -0.375d0)
H1 = 9.8298d+5*cgs_kapes*mbh**1.0d0*mdot*f(r)
rho2 = 21.921d0*alpha**(-0.7d0)*cgs_kapes**(-0.3d0)*mbh**(-0.7d0)*mdot** &
      0.4d0*r**(-1.65d0)*f(r)**0.4d0
T2 = 6.7232d+8*alpha**(-0.2d0)*cgs_kapes**0.2d0*mbh**(-0.2d0)*mdot** &
      0.4d0*r**(-0.9d0)*f(r)**0.4d0
H2 = 4639.5d0*alpha**(-0.1d0)*cgs_kapes**0.1d0*mbh**0.9d0*mdot**0.2d0*r &
      **1.05d0*f(r)**0.2d0
rho3 = 5.9458d+5*alpha**(-0.7d0)*kappa_ff**(-0.15d0)*mbh**(-0.7d0)*mdot &
      **0.55d0*r**(-1.875d0)*f(r)**0.55d0
T3 = 7.4475d+5*alpha**(-0.2d0)*kappa_ff**0.1d0*mbh**(-0.2d0)*mdot**0.3d0 &
      *r**(-0.75d0)*f(r)**0.3d0
H3 = 154.42d0*alpha**(-0.1d0)*kappa_ff**0.05d0*mbh**0.9d0*mdot**0.15d0*r &
      **1.125d0*f(r)**0.15d0
r12 = 164.17d0*alpha**0.095238d0*cgs_kapes**0.85714d0*mbh**0.095238d0* &
      mdot**0.7619d0*f(r12)**0.7619d0
r23 = 5.0556d+19*cgs_kapes**1.3333d0*kappa_ff**(-0.66667d0)*mdot** &
      0.66667d0*f(r23)**0.66667d0
