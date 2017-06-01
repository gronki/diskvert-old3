module globals

    use slf_cgs
    use slf_rk4integr
    use slf_kramers
    use slf_space
    use iso_fortran_env
    use iso_c_binding

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
    real (real64) :: max_iteration_error = 1e-6

    !   gora przedzialu obliczeniowego
    real (real64) :: htop = 60

   !   czy rozwiazywac wolniejsza metoda dajaca wszystkie rozwiazania
   logical :: cfg_balance_multi = .false.


    real(real64) :: m_bh
    real(real64) :: acc_rate
    real(real64) :: r_calc
    real(real64) :: abun_X = 0.7, abun_Z = 0.02

    real(real64) :: kram_es
    real(real64) :: kram_ff
    real(real64) :: kram_bf
    real(real64) :: kram_abs
    real(real64) :: radius
    real(real64) :: rgraw, rschw
    real(real64) :: omega
    real(real64) :: acc_edd
    real(real64) :: flux_acc
    real(real64) :: temp_eff
    real(real64) :: zscale
    real(real64) :: csound

    real(real64), parameter :: miu = 1/2d0
    real(real64), parameter :: eta = 1/12d0

    real(real64), parameter :: pi = 4*atan(real(1,real64))


contains

    subroutine init_disk(m_bh_in,acc_rate_in,r_calc_in) bind(C)
        real(real64), intent(in), value :: m_bh_in,acc_rate_in,r_calc_in
        m_bh = m_bh_in
        acc_rate = acc_rate_in
        r_calc = r_calc_in
    end subroutine

    subroutine init_abun(abun_X_in,abun_Z_in) bind(C)
        real(real64), intent(in), value :: abun_X_in,abun_Z_in
        abun_X = abun_X_in
        abun_Z = abun_Z_in
    end subroutine

    subroutine eval_globals() bind(C)
        call eval_opacity_globals
        call eval_disk_globals
    end subroutine

    subroutine eval_opacity_globals

        kram_es = (1 + abun_X) / 5
        kram_ff = 3.68d22 * (1 - abun_Z) * (1+abun_X)
        kram_bf = 4.34d25 * abun_Z * (1+abun_X)
        kram_abs = 0

        if (cfg_opacity_bf) kram_abs = kram_abs + kram_bf
        if (cfg_opacity_ff) kram_abs = kram_abs + kram_ff

    end subroutine

    subroutine eval_disk_globals

        rgraw = cgs_graw * m_bh * cgs_msun / cgs_c**2
        rschw = 2 * rgraw
        omega = sqrt( cgs_graw * m_bh * cgs_msun / (rschw*r_calc)**3 )
        acc_edd = 4*pi*cgs_graw*(m_bh * cgs_msun) / ( cgs_c * eta * cgs_kapes )
        flux_acc = 3 * cgs_graw**2 * (m_bh * cgs_msun)**2 * acc_rate      &
            / ( 2 * cgs_c * cgs_kapes * eta * (rschw*r_calc)**3. )     &
            * ( 1 - sqrt(3 / r_calc) )
        temp_eff = ( flux_acc / cgs_stef ) ** 0.25_real64
        csound = sqrt( cgs_boltz * temp_eff / (miu * cgs_mhydr) )
        zscale = csound / omega

    end subroutine



    elemental real(real64)  function kappa_abs(rho,T)
        real(real64), intent(in) :: rho,T
        kappa_abs = kram_abs * rho / T**3.5_real64
    end function

    elemental real(real64) function kappa_abs_drho(rho,T)
        real(real64), intent(in) :: rho,T
        kappa_abs_drho = kram_abs / T**3.5_real64
    end function
    elemental real(real64) function kappa_abs_1(rho,T)
        real(real64), intent(in) :: rho,T
        kappa_abs_1 = kram_abs / T**3.5_real64
    end function

    elemental real(real64) function kappa_abs_dT(rho,T)
        real(real64), intent(in) :: rho,T
        kappa_abs_dT = - 3.5_real64 * kram_abs * rho / T**(4.5_real64)
    end function
    elemental real(real64) function kappa_abs_2(rho,T)
        real(real64), intent(in) :: rho,T
        kappa_abs_2 = - 3.5_real64 * kram_abs * rho / T**(4.5_real64)
    end function

    elemental function kappa_cond(rho,T) result(kp)
        real(real64), intent(in) :: rho, T
        real(real64) :: kp
        kp = 0
    end function
    elemental function kappa_cond_1(rho,T) result(kp)
        real(real64), intent(in) :: rho, T
        real(real64) :: kp
        kp = 0
    end function
    elemental function kappa_cond_2(rho,T) result(kp)
        real(real64), intent(in) :: rho, T
        real(real64) :: kp
        kp = 0
    end function

    subroutine grid(type,h,z,n) bind(C)
        real(real64), intent(in), value :: h
        integer, intent(in), value :: n
        integer, intent(in), value :: type
        real(real64), dimension(n), intent(out) :: z
        integer :: i,aa,ab
        real(real64) :: x
        interface f_type
            pure real(real64) function f_type(x,h)
                import real64
                real(real64), intent(in) :: x,h
            end function
        end interface
        procedure(f_type), pointer :: f

        aa = 1
        ab = -1
        if ( iand(type,GRID_REVERSE_MASK) /= 0 ) then
            ! reverse grid
            aa = -1
            ab = n
        end if

        select case (iand(type,GRID_TYPE_MASK))
        case (GRID_LINEAR)
            f => f_linear
        case (GRID_LOG)
            f => f_log
        case (GRID_ASINH)
            f => f_asinh
        case default
            error stop "Incorrect grid type!"
        end select

        do i=1,n
            x = (i * aa + ab) / real(n - 1, real64)
            z(i) = f(x,h) * zscale
        end do

    contains

        pure real(real64) function f_linear(x,h)
            real(real64), intent(in) :: x,h
            f_linear = x * h
        end function
        pure real(real64) function f_log(x,h)
            real(real64), intent(in) :: x,h
            f_log = (1 + h) ** x - 1
        end function
        pure real(real64) function f_asinh(x,h)
            real(real64), intent(in) :: x,h
            f_asinh = sinh( x * asinh(h) )
        end function

    end subroutine


end module
