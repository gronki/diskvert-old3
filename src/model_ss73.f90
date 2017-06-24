module model_ss73

    use iso_fortran_env
    use iso_c_binding

    use precision
    use globals
    use slf_cgs
    use slf_space
    use slf_rk4integr
    use slf_eulerintegr
    use slf_threshold, only: thrtanh

    integer, parameter :: ny = 4,   &
        &   c_Pgas = 1, c_Prad = 2, c_Frad = 3, c_tau = 4
    integer, parameter :: na = 8,   &
        &   c_rho = 3, &
        &   c_Tgas = 2, &
        &   c_Trad = 4, &
        &   c_zscal = 1, &
        &   c_heat = 5, &
        &   c_heat_max = 6, &
        &   c_compW = 7, &
        &   c_compY = 8

contains

    subroutine init_ss73(alpha_in) bind(C)
        real(fp), intent(in), value :: alpha_in
        alpha = alpha_in
    end subroutine

    subroutine run_ss73(z,nz,y,dy,a,nmax) bind(C)
        ! ilosc przedzialow obliczeniowych
        integer(c_int), intent(in), value :: nz
        real(fp), intent(inout), dimension(nz) :: z
        real(fp), intent(inout), dimension(ny,nz) :: y,dy
        real(fp), intent(inout), dimension(na,nz) :: a
        integer(c_int), intent(out) :: nmax
        real(fp) :: h_hi, h_lo, h
        integer :: i,it


        h_hi = 12d0 * 16
        h_lo = 12d0 / 16

        write (*,'(A10,A10,A12,A7,A5)') "ITER", "H", "FLX", "MAX", "ACTN"

        do it=1,56
            h = sqrt(h_hi*h_lo)

            call grid(ior(cfg_grid,GRID_REVERSE),h,z,nz)

            y(c_Frad,1) = flux_acc
            y(c_Prad,1) = flux_acc * 2 / (3 * cgs_c)
            y(c_Pgas,1) = temp_eff * 1e-16 * cgs_k_over_mh
            y(c_tau,1)  = epsilon(1d0)

            if ( cfg_euler_integration ) then
                call eulerintegr(z, y(:,1), f, na, y, dy, a, nmax)
            else
                call rk4integr(z, y(:,1), f, na, y, dy, a, nmax)
            end if

            if ( nmax < nz ) then
                h_hi = h
                write (*,'(I10,F10.3,Es12.4,F7.1,A5)') it, h, y(c_Frad,nmax) / flux_acc, 100.0 * nmax / nz, '\/'
            else if (y(c_Frad,nz) > 1e-10 * flux_acc) then
                h_lo = h
                write (*,'(I10,F10.3,Es12.4,F7.1,A5)') it, h, y(c_Frad,nmax) / flux_acc, 100.0 * nmax / nz, '/\'
            else
                write (*,'("iteration finished")')
                exit
            end if
        end do
    end subroutine


    subroutine f(z,y,dy,a,abort)

        real(fp), intent(in) :: z, y(:)
        real(fp), intent(inout) :: dy(size(y)), a(:)
        logical, intent(inout) :: abort

        real(fp) :: compsw
        real(fp), parameter :: toler = 1
        real(fp) :: kabs, epsi, taues

        if ( y(c_Frad) < 0 ) then
            abort = .TRUE.
        end if

        a(c_zscal)  = z / zscale
        a(c_Trad)   = ( (3 * cgs_c) / (4 * cgs_stef) * y(c_Prad)) ** 0.25d0
        a(c_heat)   = alpha * omega * ( y(c_Pgas) + y(c_Prad) )


        a(c_heat_max) = 12 * cgs_kapes * y(c_Prad)  &
                & * ( y(c_Pgas) * miu * cgs_mhydr ) &
                & / ( cgs_mel * cgs_c )

        a(c_compW) = a(c_heat) / ( a(c_heat_max) - a(c_heat) )
        if ( a(c_compW) < 0 )   a(c_compW) = 0

        select case (cfg_temperature_method)
        case (EQUATION_EQUILIBR)
            a(c_Tgas)   = a(c_Trad)
        case (EQUATION_COMPTON)
            a(c_Tgas) = a(c_Trad) * (1 + a(c_compW))
            a(c_rho) = y(c_Pgas) * miu / ( cgs_k_over_mh * a(c_Tgas) )

            kabs = kram_abs * a(c_rho) * a(c_Tgas)**(-3.5d0)
            epsi = kabs / (kabs+kram_es)
            taues = (1. - epsi)/sqrt(epsi)
            a(c_compY) = 4 * cgs_boltz * a(c_tgas) &
            & / ( cgs_mel * cgs_c**2 )  &
            & * max( taues, taues**2 )
            compsw = thrtanh( (log(a(c_compY)) - log(1.0))/log(1+toler) )

            a(c_Tgas) = a(c_trad) * ( 1 + compsw * a(c_compW) )
        case (EQUATION_BALANCE)
            stop "This is not yet implemented!"
        case default
            error stop "incorrect value of cfg_temperature_method"
        end select

        a(c_rho) = y(c_Pgas) * miu / ( cgs_k_over_mh * a(c_Tgas) )

        kabs = kram_abs * a(c_rho) * a(c_Tgas)**(-3.5d0)
        epsi = kabs / (kabs+kram_es)
        taues = (1. - epsi)/sqrt(epsi)
        a(c_compY) = 4 * cgs_boltz * a(c_tgas) &
            & / ( cgs_mel * cgs_c**2 )  &
            & * max( taues, taues**2 )

        dy(c_Frad)  = a(c_heat)
        dy(c_Prad)  = - cgs_kapes * a(c_rho) / cgs_c * y(c_Frad)
        dy(c_Pgas)  = - omega**2 * a(c_rho) * z - dy(c_Prad)
        dy(c_tau)   = - cgs_kapes * a(c_rho)

    end subroutine

end module
