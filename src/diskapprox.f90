module diskapprox

    use precision
    use slf_cgs
    use globals

    implicit none

    ! f is for internal use only
    private :: f

contains

!------------------------------------- F --------------------------------------!
!                       computes the 1-sqrt(3/r) factor                        !
!----------------------------------- INPUTS -----------------------------------!
!                          r: radius (dimensionless)                           !
!---------------------------------- RETURNS -----------------------------------!
!                      f-factor for r > 3, zero otherwise                      !
!------------------------------------------------------------------------------!

    elemental function f(r)
        real(fp), intent(in) :: r
        real(fp) :: f
        if ( r > 3 ) then
            f = 1 - sqrt(3/r)
        else
            f = 0
        end if
    end function

!------------------------------ DVAPX_ZONEBOUNDS ------------------------------!
!   computes the boundaries between Prad+elsct, Pgas+elsc, Pgas+kapff zones    !
!----------------------------------- INPUTS -----------------------------------!
!                  n: number of iterations (recommended: 24)                   !
!---------------------------------- OUTPUTS -----------------------------------!
!                  r12: boundary radius between zones 1 and 2                  !
!                  r23: boundary radius between zones 2 and 3                  !
!------------------------------------------------------------------------------!

    elemental subroutine dvapx_zonebounds(r12,r23,n)
        real(fp), intent(out) :: r12,r23
        integer, intent(in) :: n
        integer :: i
        r12 = 21.844d0 * alpha**0.095238d0 * kram_es**0.85714d0 &
                * m_bh**0.095238d0 * m_dot**0.7619d0
        r23 = 8.655d+18 * kram_abs**(-0.66667d0) * kram_es**1.3333d0 &
                * m_dot**0.66667d0
        do i = 1,n
            r12 = 21.844d0 * alpha**0.095238d0 * kram_es**0.85714d0 &
                    * m_bh**0.095238d0 * m_dot**0.7619d0*f(r12)**0.7619d0
            r23 = 8.655d+18 * kram_abs**(-0.66667d0) * kram_es**1.3333d0 &
                    * m_dot**0.66667d0 * f(r23)**0.66667d0
        end do
    end subroutine

!-------------------------------- DVAPX_ZONE1 ---------------------------------!
!     calculates rho, T, H for radiation pressure and electron scattering      !
!                              dominated zone (1)                              !
!----------------------------------- INPUTS -----------------------------------!
!                                  r: radius                                   !
!---------------------------------- OUTPUTS -----------------------------------!
!                           rho: central density                               !
!                             T: central temperature                           !
!                             H: disk height                                   !
!------------------------------------------------------------------------------!

    elemental subroutine dvapx_zone1(r,rho,T,H)
        real(fp), intent(in) :: r
        real(fp), intent(out) :: rho,T,H
        rho = 0.00045904d0*1.0/m_bh*r**1.5d0/(alpha*kram_es**3*m_dot**2*f(r)**2)
        T = 4.6181d+7*alpha**(-0.25d0)*kram_es**(-0.25d0)*m_bh**(-0.25d0)*r**( &
        -0.375d0)
        H = 69679.0d0*kram_es*m_bh**1.0d0*m_dot*f(r)
    end subroutine

!-------------------------------- DVAPX_ZONE2 ---------------------------------!
! calculates rho, T, H for gas pressure and electron scattering dominated zone !
!                                     (2)                                      !
!----------------------------------- INPUTS -----------------------------------!
!                                  r: radius                                   !
!---------------------------------- OUTPUTS -----------------------------------!
!                           rho: central density                               !
!                             T: central temperature                           !
!                             H: disk height                                   !
!------------------------------------------------------------------------------!

    elemental subroutine dvapx_zone2(r,rho,T,H)
        real(fp), intent(in) :: r
        real(fp), intent(out) :: rho,T,H
        rho = 7.5985d0*alpha**(-0.7d0)*kram_es**(-0.3d0)*m_bh**(-0.7d0)*m_dot** &
        0.4d0*r**(-1.65d0)*f(r)**0.4d0
        T = 2.3314d+8*alpha**(-0.2d0)*kram_es**0.2d0*m_bh**(-0.2d0)*m_dot**0.4d0* &
        r**(-0.9d0)*f(r)**0.4d0
        H = 2734.1d0*alpha**(-0.1d0)*kram_es**0.1d0*m_bh**0.9d0*m_dot**0.2d0*r** &
        1.05d0*f(r)**0.2d0
    end subroutine

!-------------------------------- DVAXP_ZONE3 ---------------------------------!
!  calculates rho, T, H for gas pressure and free-free opacity dominated zone  !
!                                     (3)                                      !
!----------------------------------- INPUTS -----------------------------------!
!                                  r: radius                                   !
!---------------------------------- OUTPUTS -----------------------------------!
!                           rho: central density                               !
!                             T: central temperature                           !
!                             H: disk height                                   !
!------------------------------------------------------------------------------!

    elemental subroutine dvapx_zone3(r,rho,T,H)
        real(fp), intent(in) :: r
        real(fp), intent(out) :: rho,T,H
        rho = 1.3855d+5*alpha**(-0.7d0)*kram_abs**(-0.15d0)*m_bh**(-0.7d0)*m_dot &
        **0.55d0*r**(-1.875d0)*f(r)**0.55d0
        T = 3.3653d+5*alpha**(-0.2d0)*kram_abs**0.1d0*m_bh**(-0.2d0)*m_dot**0.3d0 &
        *r**(-0.75d0)*f(r)**0.3d0
        H = 103.88d0*alpha**(-0.1d0)*kram_abs**0.05d0*m_bh**0.9d0*m_dot**0.15d0*r &
        **1.125d0*f(r)**0.15d0
    end subroutine

!--------------------------------- DVAPX_SEL ----------------------------------!
!     selects the zone based on given radius and computes appropriate disk     !
!                                  parameters                                  !
!----------------------------------- INPUTS -----------------------------------!
!                           r: radius                                          !
!                         r12: radius of 1-2 boundary                          !
!                         r23: radius of 2-3 boundary                          !
!---------------------------------- OUTPUTS -----------------------------------!
!                           rho: central density                               !
!                             T: central temperature                           !
!                             H: disk height                                   !
!------------------------------------------------------------------------------!

    elemental subroutine dvapx_sel(r, r12, r23, rho, T, H)
        real(fp), intent(in) :: r, r12, r23
        real(fp), intent(out) :: rho, T, H

        if (r < r12) then
            call dvapx_zone1(r, rho, T, H)
        else if (r > r23) then
            call dvapx_zone3(r, rho, T, H)
        else
            call dvapx_zone2(r, rho, T, H)
        end if
    end subroutine

!----------------------------------- DVAPX1 -----------------------------------!
! Computes disk density, temperature and height based on 3-zone approximation. !
!        Calls DVAPX_ZONEBOUNDS to calculate boundaries between zones.         !
!----------------------------------- INPUTS -----------------------------------!
!                                  r: radius                                   !
!---------------------------------- OUTPUTS -----------------------------------!
!                           rho: central density                               !
!                             T: central temperature                           !
!                             H: disk height                                   !
!------------------------------------------------------------------------------!

    pure subroutine dvapx1(r, rho, T, H)
        real(fp), intent(in) :: r
        real(fp), intent(out) :: rho, T, H
        real(fp) :: r12, r23

        call dvapx_zonebounds(r12, r23, 24)
        call dvapx_sel(r, r12, r23, rho, T, H)

    end subroutine

end module
