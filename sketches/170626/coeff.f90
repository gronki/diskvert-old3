A(1, 1) = -4*H**3*alpha*sqrt(cgs_graw)*pi*r**(-1.5d0)*rho*sqrt(sol_mass) &
      /(m_bh*sol_rschw**(3.0d0/2.0d0)) + m_bh*m_dot*sol_mdot_edd*( &
      -1.73205080756888d0*r**(-0.5d0) + 1)
A(2, 1) = -9.0d0/16.0d0*H**4*alpha*cgs_graw**(3.0d0/2.0d0)*r**(-4.5d0)* &
      rho**2*sol_mass**(3.0d0/2.0d0)*(kram_es + kram_abs*rho/T**(7.0d0/ &
      2.0d0))/(m_bh**3*sol_rschw**(9.0d0/2.0d0)) + 4*T**4*cgs_stef
A(3, 1) = H**2*cgs_graw*r**(-3.0d0)*rho*sol_mass/(m_bh**2*sol_rschw**3) &
      - 4.0d0/3.0d0*T**4*cgs_stef/cgs_c - 2*T*cgs_boltz*rho/cgs_mhydr
M(1, 1) = 4*H**3*alpha*sqrt(cgs_graw)*pi*r**(-1.5d0)*sqrt(sol_mass)/( &
      m_bh*sol_rschw**(3.0d0/2.0d0))
M(2, 1) = (9.0d0/8.0d0)*H**4*alpha*cgs_graw**(3.0d0/2.0d0)*r**(-4.5d0)* &
      rho*sol_mass**(3.0d0/2.0d0)*(kram_es + kram_abs*rho/T**(7.0d0/ &
      2.0d0))/(m_bh**3*sol_rschw**(9.0d0/2.0d0)) + (9.0d0/16.0d0)*H**4* &
      alpha*cgs_graw**(3.0d0/2.0d0)*kram_abs*r**(-4.5d0)*rho**2* &
      sol_mass**(3.0d0/2.0d0)/(T**(7.0d0/2.0d0)*m_bh**3*sol_rschw**( &
      9.0d0/2.0d0))
M(3, 1) = -H**2*cgs_graw*r**(-3.0d0)*sol_mass/(m_bh**2*sol_rschw**3) + 2 &
      *T*cgs_boltz/cgs_mhydr
M(1, 2) = 0
M(2, 2) = -63.0d0/32.0d0*H**4*alpha*cgs_graw**(3.0d0/2.0d0)*kram_abs*r** &
      (-4.5d0)*rho**3*sol_mass**(3.0d0/2.0d0)/(T**(9.0d0/2.0d0)*m_bh**3 &
      *sol_rschw**(9.0d0/2.0d0)) - 16*T**3*cgs_stef
M(3, 2) = (16.0d0/3.0d0)*T**3*cgs_stef/cgs_c + 2*cgs_boltz*rho/cgs_mhydr
M(1, 3) = 12*H**2*alpha*sqrt(cgs_graw)*pi*r**(-1.5d0)*rho*sqrt(sol_mass) &
      /(m_bh*sol_rschw**(3.0d0/2.0d0))
M(2, 3) = (9.0d0/4.0d0)*H**3*alpha*cgs_graw**(3.0d0/2.0d0)*r**(-4.5d0)* &
      rho**2*sol_mass**(3.0d0/2.0d0)*(kram_es + kram_abs*rho/T**(7.0d0/ &
      2.0d0))/(m_bh**3*sol_rschw**(9.0d0/2.0d0))
M(3, 3) = -2*H*cgs_graw*r**(-3.0d0)*rho*sol_mass/(m_bh**2*sol_rschw**3)
