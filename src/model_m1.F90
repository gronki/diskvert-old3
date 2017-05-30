module model_m1

    use slf_cgs
    use slf_kramers
    use slf_rk4integr
    use slf_threshold
    use slf_findzer
    use slf_eulerintegr
    use ieee_arithmetic

    use globals
    use energy_balance
    use settings
    use iso_c_binding
    use iso_fortran_env
    use alphadisk
    use kind

    implicit none

    ! stale zwiazane z polem magnetycznym
    real(fp) :: alpha
    real(fp) :: beta_0
    real(fp) :: dzeta

    ! gestosc i temperatura centralna dla trybu __DISKV_SINGLE__
    real(fp) :: rho_0_user
    real(fp) :: temp_0_user

    integer, parameter :: n_vals = 7, &
                v_taues =     1, &
                v_tauth =   2, &
                v_pgas =    3, &
                v_prad =    4, &
                v_pmag =    5, &
                v_flux =    6, &
                v_fgen =    7

    integer, parameter :: n_pars = 41, &
                p_rho       = 1, &
                p_temp      = 2, &
                p_Trad      = 3, &
                p_miu       = 4, &
                p_epsi      = 5, &
                p_kappa     = 6, &
                p_pradbond  = 7, &
                p_bflux     = 8, &
                p_valfv     = 9, &
                p_csound    = 10, &
                p_vffall    = 11, &
                p_vrise     = 12, &
                p_cool_net  = 13, &
                p_energy_bil = 14, &
                p_cool_compt = 15, &
                p_ptot      = 16, &
                p_pmri      = 17, &
                p_beta      = 18, &
                p_betarad   = 19, &
                p_mri       = 20, &
                p_taues     = 21,   &
                p_compY     = 22, &
                p_compsw    = 23, &
                p_rhorad    = 24, &
                p_heat_max  = 25, &
                p_compW     = 26, &
                p_vrise_t   = 27, &
                p_vrise_tt  = 28, &
                p_cool_dyfu = 29, &
                p_temp_compton = 30, &
                p_rho_crit = 31, &
                p_pgascmpt  = 32, &
                p_flxbond   = 33, &
                p_cool_crit = 34, &
                p_heat_dyfu = 35, &
                p_tbil = 36, & ! 36 to 40
                p_zscal = 41

contains

    subroutine init_m1(alpha_in,beta_0_in) bind(C)
        real(c_double), intent(in), value :: alpha_in, beta_0_in
        alpha = alpha_in
        beta_0 = beta_0_in
    end subroutine

    subroutine run_m1(z,nz,val,der,par,nmax) bind(C)
        ! ilosc przedzialow obliczeniowych
        integer(c_int), intent(in), value :: nz
        real(c_double), intent(inout), dimension(nz) :: z
        real(c_double), intent(inout), dimension(n_vals,nz) :: val,der
        real(c_double), intent(inout), dimension(n_pars,nz) :: par
        integer(c_int), intent(out), optional :: nmax
        real(c_double) :: rho_0, temp_0
        real(c_double) :: rho_0_estim, temp_0_estim
        real(c_double) :: temp_0_hi, temp_0_lo
        real(c_double) :: rho_0_hi, rho_0_lo
        integer :: iter, iter0, itmax
        character(len=1024) :: buf

        call quick_alpha_disk(alpha)
        rho_0_estim = alphadisk_rho(1)
        temp_0_estim = alphadisk_temp(1)
        write (error_unit, fmt_meta_ec) "rho_0_estim", rho_0_estim, "Central density (estimate)"
        write (error_unit, fmt_meta_ec) "temp_0_estim", temp_0_estim, "Central gas temperature (estimate)"

        rho_0_hi = rho_0_estim * 316
        rho_0_lo = rho_0_estim / 316
        itmax = 64

        do iter0=1,itmax

            rho_0 = sqrt(rho_0_lo*rho_0_hi)

            temp_0_hi = temp_0_estim * 100
            temp_0_lo = temp_0_estim / 100

            write (error_unit,"(2A5,4A14,A7,A10)") "i0", "i", "T0", "Frad", "Fbound", "error", "clip", "actn"

            do iter=1,itmax

                temp_0 = sqrt(temp_0_hi * temp_0_lo)

                call single_m1(z,nz,val,der,par,rho_0, temp_0, nmax)

                if (nmax .lt. 2)  &
                    & error stop "All interval bad. quitting"

                if (  nmax .lt. nz-1  ) then
                    if ( dzeta .gt. 1 ) then
                        temp_0_hi = sqrt(temp_0_hi*temp_0)
                        buf = "\/! TEMP"
                    else
                        temp_0_lo = temp_0
                        buf = "/\! TEMP"
                    end if
                else
                    if (  par(p_flxbond,nmax) .lt. val(v_flux,nmax)  ) then
                        temp_0_lo = temp_0
                        buf = "/\  TEMP"
                    else if ( par(p_flxbond,nmax) .gt. val(v_flux,nmax)) then
                        temp_0_hi = temp_0
                        buf = "\/  TEMP"
                    else
                        exit
                    end if
                end if

                write (error_unit,"(2I5,3Es14.6,F14.7,F7.1,A10)") iter0, iter, &
                    & temp_0, val(v_flux,nmax), par(p_flxbond,nmax), &
                    & par(p_flxbond,nmax) / val(v_flux,nmax) - 1, &
                    & 100.0 * nmax / nz, trim(buf)

                if ( abs(2*(par(p_flxbond,nmax) - val(v_flux,nmax) )        &
                        /(par(p_flxbond,nmax) + val(v_flux,nmax) ))     &
                        .lt. max_iteration_error ) then
                    write (error_unit,"(A,Es9.1,A)") "Maximum precision ",max_iteration_error," reached."
                    exit
                end if

            end do

            if ( val(v_fgen,nmax) .gt. flux_acc ) then ! .or. nmax .lt. nz-1
                !too much flux
                rho_0_hi = rho_0
                buf = "\/ RHO"
            else if ( val(v_fgen,nmax) .lt. flux_acc ) then
                ! too lil flux
                rho_0_lo = rho_0
                buf = "/\ RHO"
            else
                buf = "== RHO"
            end if

            write (error_unit,*) "   ========================================="
            write (error_unit,"(A5,A5,3A14,A10)") "i0", "*", "rho0", "Fgen", "Facc" , "ACTN"
            write (error_unit,"(I5,A5,3Es14.6,A10)") iter0, "*", rho_0, val(v_fgen,nmax), flux_acc , trim(buf)
            write (error_unit,*) "   ========================================="

            if ( abs(2*(val(v_fgen,nmax) - flux_acc)/(val(v_fgen,nmax) + flux_acc)) .lt. max_iteration_error ) then
                write (error_unit,"(A,Es9.1,A)") "Maximum precision ",max_iteration_error," reached."
                exit
            end if

        end do
    end subroutine

    subroutine single_m1(z,nz,val,der,par,rho_0,temp_0,nmax) bind(C)
        ! ilosc przedzialow obliczeniowych
        integer(c_int), intent(in), value :: nz
        ! wspolrzedna
        real(c_double), intent(inout), dimension(nz) :: z
        ! wyniki
        real(c_double), intent(inout), dimension(n_vals,nz) :: val,der
        real(c_double), intent(inout), dimension(n_pars,nz) :: par
        ! wartosci na rowniku
        real(c_double), intent(in), value :: rho_0, temp_0
        ! na wyjsciu: przedzial, na ktorym algorytm sie wywalil
        integer(c_int), intent(out) :: nmax

        real(fp) :: Trad0

        ! make the initial conditions
        val(v_pgas,1) = cgs_k_over_mh * rho_0 * temp_0 / miu
        val(v_pmag,1) = val(v_pgas,1) / beta_0
        Trad0 = temp_0
        ! if ( cfg_temperature_method == EQUATION_BALANCE )    then
        !     call findnearzer( Trad0, fun2, (/ rho_0, temp_0 /) )
        !     if ( .not.ieee_is_normal(Trad0) ) then
        !         write (error_unit,'(A)') "Warning: temperature solver at the" &
        !             & // "equator didn't converge. Assuming the " &
        !             & // "radiation temperature!"
        !         Trad0 = temp_0
        !     end if
        ! end if
        val(v_prad,1)   = cgs_a / 3 * Trad0**4
        val(v_flux,1)   = 0
        val(v_taues,1)  = 0
        val(v_tauth,1)  = 0
        val(v_fgen,1)   = 0
        dzeta = 0.5 * alpha * (beta_0 / (val(v_pgas,1)/val(v_prad,1)) + beta_0 + 1.)

        !solve the equations
        if ( cfg_euler_integration ) then
            call eulerintegr(z, val(:,1), fder, n_pars, val, der, par, nmax)
        else
            call rk4integr(z, val(:,1), fder, n_pars, val, der, par, nmax)
        end if

        ! shift the optical depth so we have 0 on the surface
        val(v_taues,:) = val(v_taues,:) - val(v_taues,nmax) + epsilon(1d0)
        val(v_tauth,:) = val(v_tauth,:) - val(v_tauth,nmax) + epsilon(1d0)

    end subroutine


    subroutine fder(x,y,pder,a,abort)

        real(fp), intent(in) :: x, y(:)
        real(fp), intent(inout) :: pder(size(y)), a(:)
        logical, intent(inout) :: abort

        ! this constants decides about tolerance margin to MRI enable/shutoff
        ! it allows weighting parameter to vary in percent range indicated by
        ! this constant.
        real(fp), parameter :: toler = 3

        logical :: isnormal(size(a))

        real(fp) :: kabs
        real(fp) :: epsi
        real(fp) :: heat
        real(fp) :: alpha_eff
        real(fp) :: xxa,xxb
        integer :: i,n


        a(p_miu) = 0.5
        a(p_zscal) = x / zscale

        if ( (y(v_prad) .lt. 0) &
                & .or. (y(v_flux) .lt. 0) ) then
            abort = .true.
            return
        end if

        a(p_Trad) = (3 / cgs_a * y(v_prad))**0.25
        a(p_rhorad) = y(v_pgas) * a(p_miu) * cgs_mhydr / ( cgs_boltz * a(p_Trad) )


        if (a(p_temp) .eq. 0) then
            a(p_temp) = a(p_Trad)
            a(p_rho) = a(p_rhorad)
        end if


        alpha_eff = alpha
        if (cfg_allow_mri_shutdown) alpha_eff = a(p_mri) * alpha

        a(p_vrise) = omega * dzeta * x
        a(p_vrise_t) = omega * dzeta
        a(p_vrise_tt) = 0
        a(p_vffall) = omega*x
        a(p_bflux) = 2 * y(v_pmag) * a(p_vrise)
        a(p_ptot) = y(v_pgas) + y(v_pmag) + y(v_prad)
        A(p_flxbond) = 3 * cgs_c * y(v_prad) / 2
        a(p_pradbond) = 2 * y(v_flux) / ( 3 * cgs_c )

        a(p_temp) = a(p_Trad)

        heat = 2 * y(v_pmag) * a(p_vrise_t) - alpha_eff * omega * a(p_ptot)

        call calc_compW(heat, y(v_pgas), y(v_prad), a(p_heat_max), a(p_compW))
        a(p_temp_compton) = heat * cgs_mel * cgs_c &
            / ( 12 * kram_es * cgs_boltz * y(v_prad) * a(p_rho) )

        select case (cfg_temperature_method)
        case (EQUATION_EQUILIBR)
        case (EQUATION_COMPTON)
            a(p_temp) = a(p_Trad) * (1 + a(p_compW))
            call calc_compswitch(y(v_pgas),a(p_temp),a(p_compsw))
            a(p_temp) = a(p_Trad) * (1 + a(p_compsw)*a(p_compW))
        case (EQUATION_BALANCE)
            fcool_heat = heat
            fcool_prad = y(v_prad)
            fcool_pgas = y(v_pgas)
            fcool_trad = a(p_trad)
            if ( .not.cfg_balance_multi ) then
                a(p_tbil) = log(a(p_trad)**2 + a(p_temp_compton)**2) / 2
                call findzer(a(p_tbil), log(1d5), log(1d10), 1d-9, fcool)
                a(p_temp) = exp(a(p_tbil))
            else
                a(p_tbil+0:p_tbil+4) = 0

                call solve_balance_multi(a(p_tbil+0:p_tbil+4), n,&
                        & log(1d6), log(1d10), 50, fcool)

                a(p_temp) = exp(a(p_tbil+0))
            end if
        case default
            error stop "incorrect value of cfg_temperature_method"
        end select

        a(p_rho) =  y(v_pgas) * a(p_miu) * cgs_mhydr / ( cgs_boltz * a(p_temp) )

        kabs = kappa_abs(a(p_rho), a(p_temp))
        a(p_kappa) = kabs + kram_es
        epsi = kabs / a(p_kappa)
        a(p_epsi) = epsi
        a(p_taues) = (1 - epsi) / sqrt(epsi)
        a(p_compY) = 4 * cgs_boltz * a(p_temp) / ( cgs_mel * cgs_c**2 )     &
                * max( a(p_taues), a(p_taues)**2 )

        a(p_csound) = sqrt( y(v_pgas) / a(p_rho) )
        a(p_pmri) = 0.5 * a(p_csound) * a(p_rho) * omega * (r_calc*rschw)
        a(p_mri) = thrtanh((log(a(p_pmri)) - log(y(v_pmag)))/log(1+toler))

        a(p_cool_net) = kabs * (4*cgs_stef*a(p_temp)**4 - 3*y(v_prad)*cgs_c)
        if ( cfg_compton_term ) then
            a(p_cool_net) = a(p_cool_net)  &
                    & + 12 * kram_es * y(v_prad) * cgs_k_over_mec  &
                    & * (a(p_temp) - a(p_trad))
        end if
        a(p_cool_net) = a(p_cool_net) * a(p_rho)
        a(p_cool_compt) = 12 * a(p_rho) * kram_es * y(v_prad) * cgs_k_over_mec &
            * (a(p_temp) - a(p_trad))
        a(p_cool_dyfu) = kabs * a(p_rho) * 4*cgs_stef*a(p_temp)**4
        a(p_heat_dyfu) = kabs * a(p_rho) * 3*y(v_prad)*cgs_c

        if ( heat .ne. 0 ) then
            a(p_energy_bil) = 1 - a(p_cool_net) / heat
        end if

        xxa = 16 * cgs_stef * a(p_temp)**4 * kappa_abs(a(p_rho),a(p_temp))  &
            + (4*cgs_stef*a(p_temp)**4 - 3*y(v_prad)*cgs_c) &
            * ( a(p_temp) * kappa_abs_dT(a(p_rho),a(p_temp)) - a(p_rho) * kappa_abs_drho(a(p_rho),a(p_temp)) ) &
            + 12 * kram_es * y(v_prad) * cgs_k_over_mec  * a(p_temp)
        a(p_cool_crit) = xxa * a(p_rho)

        a(p_beta) = y(v_pgas) / y(v_pmag)
        a(p_betarad) =  y(v_pgas) / y(v_prad)

        a(p_valfv) = sqrt( 2 * y(v_pmag) / a(p_rho) )

        isnormal = ieee_is_normal(a)
        isnormal(p_rho) = isnormal(p_rho) .and. ( a(p_rho) > 0 )
        isnormal(p_temp) = isnormal(p_temp) .and. ( a(p_temp) > 0 )

        if ( any(.not. isnormal) ) then
            abort = .true.
            iterate_search_error: do i = 1, size(a)
                if (.not.isnormal(i)) then
                    write (0, "(Es11.3,' found in field ',I0)") a(i),i
                end if
            end do iterate_search_error
        end if


        ! alpha may be deleted by MRI shutoff
        alpha_eff = alpha
        if (cfg_allow_mri_shutdown) alpha_eff = a(p_mri) * alpha

        ! optical path
        pder(v_taues) = - kram_es * a(p_rho)
        pder(v_tauth) = - a(p_kappa) * a(p_rho) * sqrt(a(p_epsi))

        ! magnetic energy dissipation
        pder(v_flux) = 2 * y(v_pmag) * a(p_vrise_t) - alpha_eff * omega * a(p_ptot)
        if ( pder(v_flux) < 0 ) pder(v_flux) = 0

        ! magnetic pressure gradient
        if ( a(p_vrise) > 0 ) then
            pder(v_pmag) = - pder(v_flux) / a(p_vrise)
        else
            pder(v_pmag) = 0
        end if

        ! total generated flux (just for control)
        pder(v_fgen) =  alpha_eff * omega * a(p_ptot)

        ! radiative pressure gradient
        pder(v_prad) = -a(p_kappa)*a(p_rho) * (y(v_flux)  / cgs_c)

        !density gradient
        pder(v_pgas) = - (omega**2 * a(p_rho) * x + pder(v_prad) + pder(v_pmag))

    end subroutine

    subroutine solve_balance_multi(x, nx, xlo, xhi, nb, f)

        real(fp), intent(inout) :: x(:)
        real(fp), intent(in) :: xlo,xhi
        integer, intent(in) :: nb
        integer, intent(out) :: nx

        interface
            pure subroutine f(x,y,dy)
                import fp
                real(fp), intent(in) :: x
                real(fp), intent(out) :: y
                real(fp), intent(out), optional :: dy
            end subroutine
        end interface

        real(fp) :: y(nb), dy(nb), x0(nb)
        integer :: i
        real(fp) :: delx,yx,ym,xm


        do concurrent (i=1:nb)
            x0(i) = (i - 1d0)/(nb - 1) * (xhi - xlo) + xlo
        end do

        call f(x0(1),y(1))
        nx = 0
        scan_for_solutions: do i = 2, nb
            call f(x0(i),y(i))
            if ( y(i-1) * y(i) < 0 ) then
                xm = 0.5 * (x0(i-1) + x0(i))
                if ( fcool_heat /= 0 ) then
                    call f(xm,ym,yx)
                    delx = fcool_heat / yx * 1e-6
                else
                    delx = xm * 1e-9
                end if
                call linsect(x(nx+1), x0(i-1), x0(i), delx)
                nx = nx + 1
            end if
            if ( nx >= size(x) ) exit scan_for_solutions
        end do scan_for_solutions

    contains

        subroutine linsect(x,xlo0,xhi0,delx)
            real(fp), intent(in) :: xhi0,xlo0,delx
            real(fp), intent(out) :: x
            real(fp) :: xhi,xlo,yhi,ylo,y,s
            integer :: i,n

            n = int(log(abs((xhi0-xlo0)/delx))/log(2d0))

            xlo = xlo0
            xhi = xhi0

            sect_loop: do i = 1, n
                call f(xlo,ylo)
                call f(xhi,yhi)
                s = sign(1d0,yhi-ylo)

                x = ( yhi * xlo - ylo * xhi ) / ( yhi - ylo )
                call f(x,y)

                if ( y * s > 0 ) then
                    xhi = x
                elseif ( y * s < 0 ) then
                    xlo = x
                else
                    exit sect_loop
                end if

                if ( abs(xhi-xlo) < abs(delx) ) exit sect_loop
            end do sect_loop

        end subroutine

    end subroutine


    pure subroutine fun2(x,y,a)
        real(fp), intent(in) :: x,a(:)
        real(fp), intent(out) :: y
        real(fp) :: T0, rho0, Trad, kabs
        Trad = x
        rho0 = a(1)
        T0 = a(2)
        kabs = kappa_abs(rho0,T0)
        y = kabs * ( T0**4 - Trad**4 )
        if ( cfg_compton_term ) then
            y = y + kram_es * Trad**4 * cgs_k_over_mec2 * 4 * ( T0 - Trad )
        end if
    end subroutine

end module
