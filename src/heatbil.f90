module heatbalance


  use slf_cgs
  use globals, only: kappabs, kappsct, miu
  use iso_fortran_env, only: r64 => real64
  implicit none

contains

  pure subroutine heatbil(T, Trad, Pgas, heat)

    real(r64), intent(inout) :: T
    real(r64), intent(in) :: Trad, Pgas, heat
    real(r64), dimension(3,2) :: kap
    real(r64), dimension(2,2) :: F
    real(r64) :: rho, bil, bil_T, rhotemp
    integer :: i
    integer, parameter :: niter = 64

    rhotemp = Pgas * miu * cgs_mhydr / cgs_boltz

    converge: do i = 1,niter

      rho = rhotemp / T

      call kappabs(rho, T, kap(1,1), kap(2,1), kap(3,1))
      call kappsct(rho, T, kap(1,2), kap(2,2), kap(3,2))

      F(1,:) = kap(1,:)
      F(2,:) = kap(3,:) - rho / T * kap(2,:)

      bil = (-T*heat + 4*cgs_stef*rhotemp*(T - Trad)*(4*Trad**4* &
      cgs_k_over_mec2*F(1, 2) + (T + Trad)*(T**2 + Trad**2)*F(1, 1)))/( &
      T*Trad**4)
      bil_T = 4*cgs_stef*rhotemp*(T**5*F(2, 1) + 3*T**4*F(1, 1) + 4*T**2*Trad &
            **4*cgs_k_over_mec2*F(2, 2) - 4*T*Trad**5*cgs_k_over_mec2*F(2, 2 &
            ) - T*Trad**4*F(2, 1) + 4*Trad**5*cgs_k_over_mec2*F(1, 2) + Trad &
            **4*F(1, 1))/(T**2*Trad**4)


      T = T - ramp(i,niter) * bil / bil_T

    end do converge

  contains

    elemental function ramp(i,n) result(y)
      integer, intent(in) :: i,n
      real(r64) :: t,y
      t = merge(real(i) / n, 1.0, i .le. n)
      y = (3 - 2*t) * t**2
    end function

  end subroutine

end module
