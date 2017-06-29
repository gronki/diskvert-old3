program dv_mag

    use iso_fortran_env, only: r64 => real64
    use confort
    use modelmag
    use settings
    use globals
    use results
    use grid
    use fileunits
    implicit none


    character (2**10) :: buf

    integer :: errno, nmax
    integer :: i
    real(r64), dimension(:), allocatable :: z, x
    real(r64), dimension(:,:), allocatable :: val, der, par
    type(config) :: cfg

    character(len=12) :: p_labels(n_pars), v_labels(n_vals)
    real(r64) :: radius = 10, alpha = 0.1, zeta = 0.2, rho_0_user, temp_0_user

    call init_labels

    call rdargv_global
    call rdargv_rk4
    call mincf_read(cfg)
    call rdconf_global(cfg)
    call rdconf(cfg)
    call mincf_free(cfg)

    open(newunit = ulog, file = trim(outfn) // ".log", action = "write")
    open(newunit = upar, file = trim(outfn) // ".txt", action = "write")
    call init_m1(mbh, mdot, radius, alpha, zeta)

    allocate( z(ngrid) )
    forall (i = 1:ngrid)  z(i) = space_linlog(i,ngrid,htop) * fzscale(mbh, mdot, radius)

    allocate(   &
        val(n_vals,ngrid), &
        der(n_vals,ngrid), &
        par(n_pars,ngrid))

    if ( .not.cfg_single_run ) then

        call run_m1(z,val,der,par,nmax)

    else

        call single_m1(z, val, der, par, &
                & rho_0_user, temp_0_user, nmax)

    end if


    write (upar, fmparec) "alpha", alpha, "Parametr alfa"

    write (upar, fmparec) "rho_0", par(p_rho,1), "Central density"
    write (upar, fmparec) "temp_0", par(p_temp,1), "Central gas temperature"
    write (upar, fmparec) "T_rad_0", par(p_Trad,1), &
                                & "Central radiation temperature"

    write (upar, fmparec) "beta_0", val(v_pgas,1) / val(v_pmag,1), "Beta on the equator"
    write (upar, fmparec) "beta_rad_0", val(v_pgas,1) / val(v_prad,1), &
                                & "Radiative beta at the equator"
    write (upar, fmparec) "dzeta_0", zeta, "Parametr zeta"

    write (upar, fmparec) "flux_gen_top", val(v_fgen,nmax), "erg/cm2/s"
    write (upar, fmparec) "flux_rad_top", val(v_flux,nmax), "erg/cm2/s"
    write (upar, fmparec) "flux_bound_top", par(p_flxbond,nmax), "erg/cm2/s"

    close(upar)

    call write_results_3(z,ngrid,val,der,v_labels,n_vals,par,p_labels,n_pars,outfn)
    close(ulog)

contains


    subroutine rdconf(cfg)
        type(config), intent(inout) :: cfg
        integer :: errno
        character(len=2048) :: buf
        integer :: i

        call mincf_get(cfg, "zeta", buf)
        if ( cfg % not_found() )  then
            error stop "Zeta parameter between 0 and 1 " &
                & // "(key: zeta) is REQUIRED!"
        end if
        read (buf,*) zeta

        call mincf_get(cfg, "alpha", buf)
        if ( cfg % not_found() )  then
            error stop "Magnetic alpha-parameter (key: alpha) is REQUIRED!"
        end if
        read (buf,*) alpha

        call mincf_get(cfg, "radius", buf)
        if ( cfg % not_found() )  then
            error stop "Radius (key: radius) is REQUIRED!"
        end if
        read (buf,*) radius

        if ( cfg_single_run ) then

            call mincf_get(cfg, "rho0", buf)
            if ( cfg % not_found() )  then
                error stop "Midplane density (key: rho0) is REQUIRED!"
            end if
            read (buf,*) rho_0_user

            call mincf_get(cfg, "temp0", buf)
            if ( cfg % not_found() )  then
                error stop "Midplane temperature (key: temp0) is REQUIRED!"
            end if
            read (buf,*) temp_0_user

        end if

    end subroutine


    subroutine init_labels
        v_labels(v_pgas ) = 'PGAS'
        v_labels(v_prad) = 'PRAD'
        v_labels(v_pmag) = 'PMAG'
        v_labels(v_flux) = 'FRAD'
        v_labels(v_taues ) = 'TAU'
        v_labels(v_fgen) = 'FGEN'
        v_labels(v_tauth) = 'TAUTH'


        p_labels(p_mri    ) = 'MRI_ON'
        p_labels(p_compY  ) = 'COMPY'
        p_labels(p_compsw ) = 'COMP_ON'
        p_labels(p_taues  ) = 'TAUSCT'
        p_labels(p_rho   ) = 'RHO'
        p_labels(p_beta   ) = 'BETAMAG'
        p_labels(p_kappa  ) = 'KAPPA'
        p_labels(p_bflux  ) = 'FMAG'
        p_labels(p_pmri   ) = 'PMRI'
        p_labels(p_betarad) = 'BETARAD'
        p_labels(p_csound ) = 'CSOUND'
        p_labels(p_valfv  ) = 'VALVF'
        p_labels(p_vffall ) = 'VFFALL'
        p_labels(p_vrise  ) = 'VRISE'
        p_labels(p_ptot   ) = 'PTOT'
        p_labels(p_miu  ) = 'MIU'
        p_labels(p_temp  ) = 'TEMP'
        p_labels(p_epsi  ) = 'EPSI'
        p_labels(p_cool_compt  ) = 'COOL_COMPT'
        p_labels(p_vrise_t  ) = 'VRISE_T'
        p_labels(p_vrise_tt  ) = 'VRISE_TT'
        p_labels(p_pgascmpt  ) = 'PGASCMPT'
        p_labels(p_Trad  ) = 'TRAD'
        p_labels(p_flxbond  ) = 'FLXBOND'
        p_labels(p_cool_dyfu  ) = 'COOL_DYFU'
        p_labels(p_heat_dyfu  ) = 'HEAT_DYFU'
        p_labels(p_cool_net  ) = 'COOL_TOT'
        p_labels(p_cool_crit  ) = 'COOL_CRIT'
        p_labels(p_rho_crit  ) = 'RHO_CRIT'
        p_labels(p_energy_bil  ) = 'ENERGBIL'
        p_labels(p_temp_compton  ) = 'TCOMPT'
        p_labels(p_pradbond  ) = 'PRADBOND'
        p_labels(p_heat_max  ) = 'HEAT_MAX'
        p_labels(p_compW  ) = 'COMPW'
        p_labels(p_rhorad  ) = 'RHORAD'
        p_labels(p_tbil+0) = 'TBIL1'
        p_labels(p_tbil+1) = 'TBIL2'
        p_labels(p_tbil+2) = 'TBIL3'
        p_labels(p_tbil+3) = 'TBIL4'
        p_labels(p_tbil+4) = 'TBIL5'
        p_labels(p_zscal) = 'h'

    end subroutine


end program
