module globals

    use iso_fortran_env
    use iso_c_binding

    use precision
    use slf_cgs
    use slf_rk4integr
    use slf_kramers
    use slf_space

    implicit none

    !   liczba przedziałóœ
    integer :: ngrid = 2**12

    !   nazwa pliku wyjściowego
    character(len=256) :: outfn_cmdline = ""
    character(len=256) :: outfn_config = ""
    character(len=256) :: outfn = "out"

    !   Rodzaj przyjmowanej siatki
    integer, parameter ::   GRID_LINEAR = ishft(1,1), &
                        &   GRID_LOG    = ishft(2,1), &
                        &   GRID_ASINH  = ishft(3,1), &
                        &   GRID_REVERSE = 1
    INTEGER, PARAMETER :: GRID_REVERSE_MASK = 1
    INTEGER, PARAMETER :: GRID_TYPE_MASK = not(GRID_REVERSE_MASK)
    integer :: cfg_grid = GRID_LINEAR

    !   którego równania użyc?
    integer, parameter ::   EQUATION_EQUILIBR   = 0, &
                        &   EQUATION_COMPTON    = 1, &
                        &   EQUATION_BALANCE    = 2
    integer :: cfg_temperature_method = EQUATION_EQUILIBR

    !   czy w rownaniach bilansu ma byc czlon comptonowski
    logical :: cfg_compton_term = .true.

    !   czy uwzglednic przewodnictwo cieplne
    logical :: cfg_conduction = .false.

    !   czy przeliczyc raz z danymi wartosciami centralnymi czy iterowac
    !   dla spelneia warunkow brzegowych?
    logical :: cfg_single_run = .false.

    !   czy umozliwic wylaczenie MRI?
    logical :: cfg_allow_mri_shutdown = .false.

    !   nieprzezroczystości free-free i bound-free
    logical :: cfg_opacity_ff = .true.
    logical :: cfg_opacity_bf = .false.

    !   czy stosowac schemat eulera zamiast rk4
    logical :: cfg_euler_integration = .false.

    !   jaki jest dopuszczalny blad iteracji
    real (fp) :: max_iteration_error = 1e-6

    !   gora przedzialu obliczeniowego
    real (fp) :: htop = 60

   !   czy rozwiazywac wolniejsza metoda dajaca wszystkie rozwiazania
   logical :: cfg_balance_multi = .false.

   ! parametry globalne
    real(fp) :: m_bh
    real(fp) :: m_dot
    real(fp) :: r_calc
    real(fp) :: abun_X = 0.7381
    real(fp) :: abun_Z = 0.0134

    ! parametry modeli
    real(fp) :: alpha, zeta

    real(fp) :: kram_es = 3.45504e-01
    real(fp) :: kram_ff = 6.31050e+22
    real(fp) :: kram_bf = 1.01081e+24
    real(fp) :: kram_abs = 6.31050e+22
    real(fp) :: radius
    real(fp) :: rschw
    real(fp) :: omega
    real(fp) :: flux_acc
    real(fp) :: temp_eff
    real(fp) :: zscale
    real(fp) :: csound

    real(fp), parameter :: miu = 1/2d0

    real(fp), parameter :: pi = 4*atan(real(1,fp))


contains

    subroutine init_disk(m_bh_in,acc_rate_in,r_calc_in)   &
          bind(C, name = 'dv_init_disk')
        real(fp), intent(in), value :: m_bh_in,acc_rate_in,r_calc_in
        m_bh = m_bh_in
        m_dot = acc_rate_in
        r_calc = r_calc_in
    end subroutine

    subroutine init_abun(abun_X_in,abun_Z_in)   &
          bind(C, name = 'dv_init_abun')
        real(fp), intent(in), value :: abun_X_in,abun_Z_in
        abun_X = abun_X_in
        abun_Z = abun_Z_in
    end subroutine

    subroutine eval_globals() bind(C, name = 'dv_eval_globals')
        call eval_opacity_globals
        call eval_disk_globals
    end subroutine

    subroutine eval_opacity_globals

        kram_es = cgs_kapes_hydrogen * (1 + abun_X) / 2
        kram_ff = 3.68d22 * (1 - abun_Z) * (1 + abun_X)
        kram_bf = 4.34d25 * abun_Z * (1 + abun_X)
        kram_abs = 0

        if (cfg_opacity_bf) kram_abs = kram_abs + kram_bf
        if (cfg_opacity_ff) kram_abs = kram_abs + kram_ff

    end subroutine

    elemental subroutine cylinder(r, rschw, omega, facc)
        real(fp), intent(in) :: r
        real(fp), intent(out) :: rschw, omega, facc
        rschw = sol_rschw * m_bh
        omega = sol_omega / sqrt( m_bh**2 * r**3 )
        facc = sol_facc_edd * (m_dot / m_bh) * (1 - sqrt(3/r)) / r**3
    end subroutine

    subroutine eval_disk_globals
        call cylinder(r_calc, rschw, omega, flux_acc)
        temp_eff = ( flux_acc / cgs_stef ) ** 0.25_fp
        csound = sqrt( cgs_boltz * temp_eff / (miu * cgs_mhydr) )
        zscale = csound / omega
    end subroutine

    elemental real(fp)  function kappa_abs(rho,T)
        real(fp), intent(in) :: rho,T
        kappa_abs = kram_abs * rho / T**3.5_fp
    end function

    elemental real(fp) function kappa_abs_drho(rho,T)
        real(fp), intent(in) :: rho,T
        kappa_abs_drho = kram_abs / T**3.5_fp
    end function
    elemental real(fp) function kappa_abs_1(rho,T)
        real(fp), intent(in) :: rho,T
        kappa_abs_1 = kram_abs / T**3.5_fp
    end function

    elemental real(fp) function kappa_abs_dT(rho,T)
        real(fp), intent(in) :: rho,T
        kappa_abs_dT = - 3.5_fp * kram_abs * rho / T**(4.5_fp)
    end function
    elemental real(fp) function kappa_abs_2(rho,T)
        real(fp), intent(in) :: rho,T
        kappa_abs_2 = - 3.5_fp * kram_abs * rho / T**(4.5_fp)
    end function

    elemental function kappa_cond(rho,T) result(kp)
        real(fp), intent(in) :: rho, T
        real(fp) :: kp
        kp = 0
    end function
    elemental function kappa_cond_1(rho,T) result(kp)
        real(fp), intent(in) :: rho, T
        real(fp) :: kp
        kp = 0
    end function
    elemental function kappa_cond_2(rho,T) result(kp)
        real(fp), intent(in) :: rho, T
        real(fp) :: kp
        kp = 0
    end function


end module
