bil = 4*Pgas*cgs_mhydr*cgs_stef*miu*(T - Trad)*(4*Trad**4*cgs_boltz*F(1, &
      2)/(cgs_c**2*cgs_mel) + (T + Trad)*(T**2 + Trad**2)*F(1, 1))/(T* &
      cgs_boltz*heat) - 1
bil_T = 4*Pgas*cgs_mhydr*cgs_stef*miu*(T**5*cgs_c**2*cgs_mel*F(2, 1) + 3 &
      *T**4*cgs_c**2*cgs_mel*F(1, 1) + 4*T**2*Trad**4*cgs_boltz*F(2, 2 &
      ) - 4*T*Trad**5*cgs_boltz*F(2, 2) - T*Trad**4*cgs_c**2*cgs_mel*F( &
      2, 1) + 4*Trad**5*cgs_boltz*F(1, 2) + Trad**4*cgs_c**2*cgs_mel*F( &
      1, 1))/(T**2*cgs_boltz*cgs_c**2*cgs_mel*heat)
