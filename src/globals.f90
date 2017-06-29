module globals

    use iso_fortran_env, only: r64 => real64
    use abundance
    use slf_cgs

    implicit none !------------------------------------------------------------!

    real, parameter, private :: X0 = 0.7381
    real, parameter, private :: Z0 = 0.0134

    real(r64) :: kappa_es = cgs_kapes_hydrogen * (1 + X0) / 2
    real(r64) :: kappa_ff_0 = 3.68d22 * (1 - Z0) * (1 + X0)
    real(r64) :: kappa_bf_0 = 4.34d25 * Z0 * (1 + X0)
    real(r64) :: kappa_abs_0 = 3.68d22 * (1 - Z0) * (1 + X0)

    logical :: kramers_opacity_ff = .true.
    logical :: kramers_opacity_bf = .false.

    real(r64) :: mbh, mdot

    real(r64), parameter :: miu = 0.50
    real(r64), parameter :: pi = 4*atan(real(1,r64))

contains !---------------------------------------------------------------------!


    elemental subroutine cylinder(mbh, mdot, r, omega, flux, teff, zscale)
        real(r64), intent(in) :: mbh, mdot, r
        real(r64), intent(out) :: omega, flux, teff, zscale
        omega = sol_omega / sqrt( mbh**2 * r**3 )
        flux = sol_facc_edd * (mdot / mbh) * (1 - sqrt(3/r)) / r**3
        teff = ( flux / cgs_stef ) ** (1d0 / 4d0)
        zscale = sqrt( 2 * cgs_k_over_mh * teff ) / omega
    end subroutine

    elemental function fzscale(mbh, mdot, r) result(zscale)
        real(r64), intent(in) :: mbh, mdot, r
        real(r64) :: omega, flux, teff, zscale
        omega = sol_omega / sqrt( mbh**2 * r**3 )
        flux = sol_facc_edd * (mdot / mbh) * (1 - sqrt(3/r)) / r**3
        teff = ( flux / cgs_stef ) ** (1d0 / 4d0)
        zscale = sqrt( 2 * cgs_k_over_mh * teff ) / omega
    end function

!------------------------------------------------------------------------------!

    elemental function fksct(rho,T) result(ksct)
        real(r64), intent(in) :: rho,T
        real(r64):: ksct
        ksct = cgs_kapes_hydrogen * (1 + abun_X) / 2
        ! red = ( 1 + 2.7e11 * rho2 / T2**2 ) * ( 1 + (T2 / 4.5e8)**0.86 )
        ! ksct = ksct0 / red
    end function

    elemental subroutine KAPPSCT(rho,T,kap,krho,kT)
        real(r64), intent(in) :: rho,T
        real(r64), intent(out) :: kap,krho,kT
        kap = cgs_kapes_hydrogen * (1 + abun_X) / 2
        krho = 0
        kT = 0
    end subroutine

!------------------------------------------------------------------------------!

    elemental function kramers(X, Z) result(kab0)
        real(r64), intent(in) :: X, Z
        real(r64) :: kab0, kff0, kbf0
        kff0 = 3.68d22 * (1 - Z) * (1 + X)
        kbf0 = 4.34d25 * Z * (1 + X)
        kab0 = merge(kff0, 0.0_r64, kramers_opacity_ff)   &
             + merge(kbf0, 0.0_r64, kramers_opacity_bf)
    end function

!------------------------------------------------------------------------------!

    elemental function fkabs(rho,T) result(kabs)
        real(r64), intent(in) :: rho,T
        real(r64) :: kabs
        kabs = kramers(abun_X,abun_Z) * rho * T ** (-7d0/2d0)
    end function
    elemental function fkabs_1(rho,T) result(kabs)
        real(r64), intent(in) :: rho,T
        real(r64) :: kabs
        kabs = kramers(abun_X,abun_Z) * T ** (-7d0/2d0)
    end function
    elemental function fkabs_2(rho,T) result(kabs)
        real(r64), intent(in) :: rho,T
        real(r64) :: kabs
        kabs = (-7d0/2d0) * kramers(abun_X,abun_Z) * rho * T ** (-9d0/2d0)
    end function

!------------------------------------------------------------------------------!

    elemental subroutine KAPPABS(rho,T,kap,krho,kT)
        real(r64), intent(in) :: rho,T
        real(r64), intent(out) :: kap,krho,kT

        kap = kramers(abun_X,abun_Z) * rho * T ** (-7d0/2d0)
        krho = kap / rho
        kT = (-7d0/2d0) * kap / T
    end subroutine

!------------------------------------------------------------------------------!

    elemental real(r64) function kappa_abs_drho(rho,T) result(kappa)
        real(r64), intent(in) :: rho,T
        kappa = kappa_abs_0 / T ** (7.0_r64 / 2.0_r64)
    end function
    elemental real(r64) function kappa_abs_1(rho,T) result(kappa)
        real(r64), intent(in) :: rho,T
        kappa = kappa_abs_0 / T ** (7.0_r64 / 2.0_r64)
    end function

    elemental real(r64) function kappa_abs_dT(rho,T) result(kappa)
        real(r64), intent(in) :: rho,T
        kappa = - (7.0_r64 / 2.0_r64) * kappa_abs_0 &
            & * rho / T ** (9.0_r64 / 2.0_r64)
    end function
    elemental real(r64) function kappa_abs_2(rho,T) result(kappa)
        real(r64), intent(in) :: rho,T
        kappa = - (7.0_r64 / 2.0_r64) * kappa_abs_0 &
            & * rho / T ** (9.0_r64 / 2.0_r64)
    end function

!------------------------------------------------------------------------------!

    elemental function fkcnd(rho,T) result(kcnd)
        real(r64), intent(in) :: rho,T
        real(r64) :: kcnd
        kcnd = 0
    end function
    elemental function fkcnd_1(rho,T) result(kcnd)
        real(r64), intent(in) :: rho,T
        real(r64) :: kcnd
        kcnd = 0
    end function
    elemental function fkcnd_2(rho,T) result(kcnd)
        real(r64), intent(in) :: rho,T
        real(r64) :: kcnd
        kcnd = 0
    end function

!------------------------------------------------------------------------------!

end module
