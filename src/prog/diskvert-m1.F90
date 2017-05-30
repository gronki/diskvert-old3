program m1

    use confort

    use model_m1

    use setup
    use settings
    use globals
    use results


    implicit none

    integer :: nmax

    character (len=1024) :: buf

    integer :: meta_lun
    integer :: errno
    integer :: i
    real(fp), dimension(:), allocatable :: z, x
    real(fp), dimension(:,:), allocatable :: val, der, par

    character(len=12) :: p_labels(n_pars), v_labels(n_vals)

    errno = 0
    call init(errno, conf_reader = disk_read_local)

    call init_labels

    open(newunit = meta_lun, file = trim(outfn) // ".txt", action = "write")

    allocate( z(ngrid) )
    call grid(cfg_grid,htop,z,ngrid)


    allocate(   &
        val(n_vals,ngrid), &
        der(n_vals,ngrid), &
        par(n_pars,ngrid))

    if ( .not.cfg_single_run ) then

        call run_m1(z,ngrid,val,der,par,nmax)

    else

        call single_m1(z, ngrid, val, der, par, &
                & rho_0_user, temp_0_user, nmax)

    end if


    write (meta_lun, fmt_meta_ec) "alpha", alpha, "Parametr alfa"

    write (meta_lun, fmt_meta_ec) "rho_0", par(p_rho,1), "Central density"
    write (meta_lun, fmt_meta_ec) "temp_0", par(p_temp,1), "Central gas temperature"
    write (meta_lun, fmt_meta_ec) "T_rad_0", par(p_Trad,1), &
                                & "Central radiation temperature"

    write (meta_lun, fmt_meta_ec) "beta_0", beta_0, "Beta on the equator"
    write (meta_lun, fmt_meta_ec) "beta_rad_0", val(v_pgas,1) / val(v_prad,1), &
                                & "Radiative beta at the equator"
    write (meta_lun, fmt_meta_ec) "dzeta_0", dzeta, "Parametr dzeta"

    write (meta_lun, fmt_meta_ec) "flux_gen_top", val(v_fgen,nmax), "erg/cm2/s"
    write (meta_lun, fmt_meta_ec) "flux_rad_top", val(v_flux,nmax), "erg/cm2/s"
    write (meta_lun, fmt_meta_ec) "flux_bound_top", par(p_flxbond,nmax), "erg/cm2/s"

    call write_globals(meta_lun)
    call write_settings(meta_lun)

    close(meta_lun)

    call write_results_3(z,ngrid,val,der,v_labels,n_vals,par,p_labels,n_pars,outfn)

contains


    subroutine disk_read_local(cfg,errno)
        integer, intent(inout) :: errno
        type(confort_c), intent(inout) :: cfg
        character(len=2048) :: buf
        integer :: i

        call mincf_get(cfg, "beta_0", buf, errno)
        if ( iand(errno,MINCF_NOT_FOUND) .ne. 0 )  then
            call mincf_free(cfg)
            error stop "Magnetic beta (Pgas/Pmag) on equatorial plane " &
                & // "(key: beta_0) is REQUIRED!"
        end if
        read (buf,*) beta_0

        call mincf_get(cfg, "alpha", buf, errno)
        if ( iand(errno,MINCF_NOT_FOUND) .ne. 0 )  then
            call mincf_free(cfg)
            error stop "Magnetic alpha-parameter (key: alpha) is REQUIRED!"
        end if
        read (buf,*) alpha

        if ( cfg_single_run ) then

            call mincf_get(cfg, "rho0", buf, errno)
            if ( iand(errno,MINCF_NOT_FOUND) .ne. 0 )  then
                call mincf_free(cfg)
                error stop "Midplane density (key: rho0) is REQUIRED!"
            end if
            read (buf,*) rho_0_user

            call mincf_get(cfg, "temp0", buf, errno)
            if ( iand(errno,MINCF_NOT_FOUND) .ne. 0 )  then
                call mincf_free(cfg)
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
