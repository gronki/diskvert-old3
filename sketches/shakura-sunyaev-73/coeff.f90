subroutine coeffs(mbh,mdot,r,alpha,rho,T,H,A,M) bind(C)
implicit none
double precision, intent(in), value :: mbh,mdot,r,alpha,rho,T,H
double precision, intent(out) :: A(3,1), M(3,3)
A(1, 1) = 902184.7d0*H**3*alpha*rho*sqrt(1/(mbh**2*r**3)) - 1.975034d+18 &
      *mbh*mdot*(-1.732051d0*sqrt(1.0/r) + 1.0d0)
A(2, 1) = 2.081514d+14*H**4*alpha*rho**2*(1/(mbh**2*r**3))**1.5d0*( &
      6.13d+22*T**(-3.5d0)*rho + 0.3379305d0) - 0.0002268204d0*T**4
A(3, 1) = -5.154318d+9*H**2*rho*(1/(mbh**2*r**3))**1.0d0 + 2.521972d-15* &
      T**4 + 1.650222d+8*T*rho
M(1, 1) = 902184.7d0*H**3*alpha*sqrt(1/(mbh**2*r**3))
M(2, 1) = 1.275968d+37*H**4*T**(-3.5d0)*alpha*rho**2*(1/(mbh**2*r**3))** &
      1.5d0 + 4.163028d+14*H**4*alpha*rho*(1/(mbh**2*r**3))**1.5d0*( &
      6.13d+22*T**(-3.5d0)*rho + 0.3379305d0)
M(3, 1) = -5.154318d+9*H**2*(1/(mbh**2*r**3))**1.0d0 + 1.650222d+8*T
M(1, 2) = 0
M(2, 2) = -4.465888d+37*H**4*T**(-4.5d0)*alpha*rho**3*(1/(mbh**2*r**3)) &
      **1.5d0 - 0.0009072817d0*T**3
M(3, 2) = 1.008789d-14*T**3 + 1.650222d+8*rho
M(1, 3) = 2706554.0d0*H**2*alpha*rho*sqrt(1/(mbh**2*r**3))
M(2, 3) = 8.326056d+14*H**3*alpha*rho**2*(1/(mbh**2*r**3))**1.5d0*( &
      6.13d+22*T**(-3.5d0)*rho + 0.3379305d0)
M(3, 3) = -1.030864d+10*H*rho*(1/(mbh**2*r**3))**1.0d0
end subroutine
