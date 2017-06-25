rho1 = 0.00045904d0*1.0/mbh*r**1.5d0/(alpha*cgs_kapes**3*mdot**2*f(r)**2 &
      )
T1 = 4.6181d+7*alpha**(-0.25d0)*cgs_kapes**(-0.25d0)*mbh**(-0.25d0)*r**( &
      -0.375d0)
H1 = 69679.0d0*cgs_kapes*mbh**1.0d0*mdot*f(r)
rho2 = 7.5985d0*alpha**(-0.7d0)*cgs_kapes**(-0.3d0)*mbh**(-0.7d0)*mdot** &
      0.4d0*r**(-1.65d0)*f(r)**0.4d0
T2 = 2.3314d+8*alpha**(-0.2d0)*cgs_kapes**0.2d0*mbh**(-0.2d0)*mdot** &
      0.4d0*r**(-0.9d0)*f(r)**0.4d0
H2 = 2734.1d0*alpha**(-0.1d0)*cgs_kapes**0.1d0*mbh**0.9d0*mdot**0.2d0*r &
      **1.05d0*f(r)**0.2d0
rho3 = 1.3855d+5*alpha**(-0.7d0)*kappa_ff**(-0.15d0)*mbh**(-0.7d0)*mdot &
      **0.55d0*r**(-1.875d0)*f(r)**0.55d0
T3 = 3.3653d+5*alpha**(-0.2d0)*kappa_ff**0.1d0*mbh**(-0.2d0)*mdot**0.3d0 &
      *r**(-0.75d0)*f(r)**0.3d0
H3 = 103.88d0*alpha**(-0.1d0)*kappa_ff**0.05d0*mbh**0.9d0*mdot**0.15d0*r &
      **1.125d0*f(r)**0.15d0
r12 = 21.844d0*alpha**0.095238d0*cgs_kapes**0.85714d0*mbh**0.095238d0* &
      mdot**0.7619d0*f(r12)**0.7619d0
r23 = 8.655d+18*cgs_kapes**1.3333d0*kappa_ff**(-0.66667d0)*mdot** &
      0.66667d0*f(r23)**0.66667d0
