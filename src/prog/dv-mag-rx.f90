program dv_mag_relax

  use confort
  use globals
  use settings
  use iso_fortran_env, only: sp => real32, dp => real64
  use fileunits
  use relaxation
  use relaxutils
  use slf_deriv, only: deriv
  use rxsettings
  use ss73solution, only: apxdisk, apx_estim, apx_refin
  use grid

  !----------------------------------------------------------------------------!
  implicit none

  type(config) :: cfg
  integer :: model, errno
  integer :: ny = 3, i, iter, nitert = 0, KL, KU
  integer, dimension(3) :: niter = [ 24, 12, 30 ]
  real(dp), allocatable, target :: x(:), x0(:), Y(:), dY(:), M(:,:), YY(:,:)
  real(dp), allocatable :: MB(:,:)
  real(dp), pointer :: yv(:,:)
  integer, dimension(:), allocatable :: ipiv
  logical, dimension(:), allocatable :: errmask
  real(dp), pointer, dimension(:) :: y_rho, y_temp, y_frad, y_pmag, &
        & y_Trad, y_fcnd
  real(dp), allocatable, dimension(:) :: tau, heat
  real(dp) :: rhoc, Tc, Hdisk, err, err0, ramp, beta_0
  character(*), parameter :: fmiter = '(I5,2X,ES9.2,2X,F5.1,"%")'
  logical :: user_ff, user_bf, converged
  integer, dimension(6) :: c_
  integer, parameter :: upar = 92

  real(dp), parameter :: typical_hdisk = 15

  integer, parameter :: ncols  = 25, &
      c_rho = 1, c_temp = 2, c_trad = 3, &
      c_pgas = 4, c_prad = 5, c_pmag = 6, &
      c_frad = 7, c_fmag = 8, c_fcnd = 9, &
      c_heat = 10, c_vrise = 11, &
      c_ksct = 12, c_kabs = 13, c_kcnd = 14, &
      c_tau = 15, c_taues = 16, c_tauth = 17, &
      c_tavg = 18, c_beta = 19, c_coldens = 20, &
      c_cool0 = 21, c_coolb = 22, c_coolc = 23, &
      c_compy = 24, c_kabp = 25
  character(8), dimension(ncols) :: labels

  labels(c_rho) = 'rho'
  labels(c_temp) = 'temp'
  labels(c_trad) = 'trad'
  labels(c_pgas) = 'pgas'
  labels(c_prad) = 'prad'
  labels(c_pmag) = 'pmag'
  labels(c_frad) = 'frad'
  labels(c_fmag) = 'fmag'
  labels(c_fcnd) = 'fcnd'
  labels(c_vrise) = 'vrise'
  labels(c_heat) = 'heat'
  labels(c_vrise) = 'vrise'
  labels(c_ksct) = 'ksct'
  labels(c_kabs) = 'kabs'
  labels(c_kabp) = 'kabp'
  labels(c_kcnd) = 'kcnd'
  labels(c_tau) = 'tau'
  labels(c_taues) = 'taues'
  labels(c_tauth) = 'tauth'
  labels(c_tavg) = 'tavg'
  labels(c_beta) = 'beta'
  labels(c_coldens) = 'coldens'
  labels(c_cool0) = 'cool0'
  labels(c_coolb) = 'coolb'
  labels(c_coolc) = 'coolc'
  labels(c_compy) = 'compy'

  !----------------------------------------------------------------------------!
  ! default values

  ngrid = -1
  htop = 150

  !----------------------------------------------------------------------------!
  ! initialize the globals, read the config etc

  call rdargvgl
  call rdargvrx
  call mincf_read(cfg)
  call rdconfgl(cfg)
  call rdconf(cfg)
  call mincf_free(cfg)


  open(upar, file = trim(outfn) // '.txt', action = 'write')

  user_ff = use_opacity_ff
  user_bf = use_opacity_bf
  converged = .true.

  !----------------------------------------------------------------------------!
  ! calculate the global parameters

  call cylinder(mbh, mdot, radius, omega, facc, teff, zscale)

  ! calculate the SS73 solution to obtain the initial profile for iteration
  call apx_estim(mbh, mdot, radius, alpha, rhoc, Tc, Hdisk)
  call apx_refin(mbh, mdot, radius, alpha, rhoc, Tc, Hdisk)

  beta_0 = 2 * zeta / alpha + nu - 1

  if (beta_0 < 0) error stop "MAGNETIC BETA cannot be negative! " &
        & // "zeta > alpha / 2!!!"

  write (upar,fmpare) 'beta_0', beta_0

  if ((tgrid == GRID_LOG .or. tgrid == GRID_ASINH) &
      .and. cfg_adjust_height_beta) &
    htop = htop * (1 + 1 / sqrt(beta_0)) * (hdisk / (3 * zscale))

  !----------------------------------------------------------------------------!
  ! if the grid number has not been set, choose the default

  if (ngrid .eq. -1) then
    select case (tgrid)
    case (grid_linear)
      ngrid = ceiling(8 * htop)
    case (grid_log, grid_asinh)
      ngrid = ceiling(300 * log(1 + htop / typical_hdisk))
    case (grid_pow2)
      ngrid = ceiling(4 * htop)
    case default
      error stop "this grid is not supported"
    end select
  end if

  !----------------------------------------------------------------------------!

  ! get the model number
  model = mrx_number( 'D', .TRUE., .FALSE. )
  call mrx_sel_nvar(model, ny)
  call mrx_sel_hash(model, C_)

  allocate( x(ngrid), x0(ngrid), Y(ny*ngrid), dY(ny*ngrid),  M(ny*ngrid,ny*ngrid), tau(ngrid), heat(ngrid), YY(ncols,ngrid) )
  allocate(errmask(ny*ngrid), ipiv(ny*ngrid))

  !----------------------------------------------------------------------------!
  ! generate the grid

  generate_grid: block

    select case (tgrid)
    case (grid_linear)
      forall (i = 1:ngrid)
        x(i)  = space_linear(i, ngrid, htop) * zscale
        x0(i) = space_linear(i, ngrid, htop) / htop
      end forall
    case (grid_log)
      forall (i = 1:ngrid)
        x(i)  = space_linlog(i, ngrid, htop / typical_hdisk) &
              * typical_hdisk * zscale
        x0(i) = space_linlog(i, ngrid, htop / typical_hdisk) &
              * typical_hdisk / htop
      end forall
    case (grid_asinh)
      forall (i = 1:ngrid)
        x(i)  = space_asinh(i, ngrid, htop / typical_hdisk) &
              * typical_hdisk * zscale
        x0(i) = space_asinh(i, ngrid, htop / typical_hdisk) &
              * typical_hdisk / htop
      end forall
    case (grid_pow2)
      forall (i = 1:ngrid)
        x(i)  = space_pow2(i, ngrid, real(100, dp)) * htop * zscale
        x0(i) = space_pow2(i, ngrid, real(100, dp))
      end forall
    case default
      error stop "this grid is not supported"
    end select

  end block generate_grid

  !----------------------------------------------------------------------------!
  ! set pointers

  y_rho   => Y(C_(1)::ny)
  y_temp  => Y(C_(2)::ny)
  y_trad  => Y(C_(3)::ny)
  y_frad  => Y(C_(4)::ny)
  y_pmag  => Y(C_(5)::ny)
  YV(1:ny,1:ngrid) => Y

  !----------------------------------------------------------------------------!
  ! generate the initial disk profile

  forall (i = 1:ngrid)
    y_frad(i) = x0(i) * facc
    y_temp(i) = (1 - x0(i)) * (Tc - 0.841 * Teff) + 0.841 * Teff
    y_rho(i) =  rhoc * (exp(-0.5*(x(i)/Hdisk)**2) + 1e-6)

    y_pmag(i) = 2 * cgs_k_over_mh * y_rho(i) * y_temp(i)   &
          & / (beta_0 * exp(- 0.5 * (x(i) / Hdisk)**2 ) + 1e-2)
  end forall

  !----------------------------------------------------------------------------!
  ! do the initial relaxation with only electron scattering opacity

  if (cfg_write_all_iters) call saveiter(0)

  use_opacity_ff = .false.
  use_opacity_bf = .false.
  err0 = 0

  relx_opacity_es : do iter = 1,niter(1)

    call mrx_matrix(model, x, Y, M, dY)
    ! call dgesv(size(M,2), 1, M, size(M,1), ipiv, dY, size(dY), errno)
    call m2band(M, MB, KL, KU)
    call dgbsv(size(MB,2), KL, KU, 1, MB, size(MB,1), ipiv, dY, size(dY), errno)

    forall (i = 1:ngrid*ny) errmask(i) = (Y(i) .ne. 0) &
        & .and. ieee_is_normal(dY(i))
    err = sqrt(sum((dY/Y)**2, errmask) / count(errmask))
    ramp = ramp5(iter, niter(1) - 6)
    write(uerr,fmiter) nitert+1, err, 100*ramp

    Y(:) = Y + dY * ramp

    nitert = nitert + 1
    if (cfg_write_all_iters) call saveiter(nitert)

    if (err < 1e-5 .and. err0 / err > 5) then
      write (uerr, '("convergence reached with error = ",ES9.2)') err
      exit relx_opacity_es
    end if

    err0 = err

  end do relx_opacity_es


  !----------------------------------------------------------------------------!
  ! do relaxation with other opacities

  use_opacity_ff = user_ff
  use_opacity_bf = user_bf

  err0 = 0

  if ( use_opacity_bf .or. use_opacity_ff ) then

    write (uerr,*) '--- ff+bf opacity is on'

    relx_opacity_full : do iter = 1,niter(2)

      call mrx_matrix(model, x, Y, M, dY)
      ! call dgesv(size(M,2), 1, M, size(M,1), ipiv, dY, size(dY), errno)
      call m2band(M, MB, KL, KU)
      call dgbsv(size(MB,2), KL, KU, 1, MB, size(MB,1), ipiv, dY, size(dY), errno)

      forall (i = 1:ngrid*ny) errmask(i) = (Y(i) .ne. 0) &
          & .and. ieee_is_normal(dY(i))
      err = sqrt(sum((dY/Y)**2, errmask) / count(errmask))
      ramp = ramp2(iter, niter(2) - 2)
      write(uerr,fmiter) nitert+1, err, 100*ramp

      Y(:) = Y + dY * ramp

      nitert = nitert + 1
      if (cfg_write_all_iters) call saveiter(nitert)

      if (err < 1e-5 .and. err0 / err > 5) then
        write (uerr, '("convergence reached with error = ",ES9.2)') err
        exit relx_opacity_full
      end if

      err0 = err

    end do relx_opacity_full


  end if

  !----------------------------------------------------------------------------!
  ! if the coorna was requested, relax the gas temperature

  if ( cfg_temperature_method .ne. EQUATION_DIFFUSION ) then

    write (uerr,*) '--- corona is on'

    err0 = 0

    call mrx_transfer(model, &
        mrx_number(cfg_temperature_method, .TRUE., use_conduction), x, Y)

    call mrx_sel_nvar(model, ny)
    call mrx_sel_hash(model, c_)

    deallocate(dY, M, errmask, ipiv)
    allocate(dY(ny*ngrid), M(ny*ngrid,ny*ngrid))
    allocate(errmask(ny*ngrid), ipiv(ny*ngrid))

    y_rho   => Y(C_(1)::ny)
    y_temp  => Y(C_(2)::ny)
    y_trad  => Y(C_(3)::ny)
    y_frad  => Y(C_(4)::ny)
    y_pmag  => Y(C_(5)::ny)
    YV(1:ny,1:ngrid) => Y

    if (cfg_temperature_method == EQUATION_BALANCE) niter(3) = niter(3) + 6

    if (use_conduction) then
      niter(3) = niter(3) + 12
      y_fcnd => Y(C_(6)::ny)
      y_fcnd(:) = 0
    end if

    call deriv(x, y_frad, heat)

    forall (i = 1:ngrid)
      y_temp(i) = heat(i) * (cgs_mel * cgs_c**2) &
            & / (16 * cgs_boltz * (cgs_kapes*y_rho(i)) &
            & * (cgs_stef*y_trad(i)**4) )
      y_temp(i) = sqrt(y_trad(i)**2 + y_temp(i)**2)
    end forall

    nitert = nitert + 1
    if (cfg_write_all_iters) call saveiter(nitert)

    !--------------------------------------------------------------------------!

    relx_corona : do iter = 1,niter(3)

      call mrx_matrix(model, x, Y, M, dY)
      ! call dgesv(size(M,2), 1, M, size(M,1), ipiv, dY, size(dY), errno)
      call m2band(M, MB, KL, KU)
      call dgbsv(size(MB,2), KL, KU, 1, MB, size(MB,1), ipiv, dY, size(dY), errno)

      if (errno .ne. 0) exit relx_corona

      forall (i = 1:ngrid*ny) errmask(i) = (Y(i) .ne. 0) &
          & .and. ieee_is_normal(dY(i))
      err = sqrt(sum((dY/Y)**2, errmask) / count(errmask))
      ramp = ramp5(iter, niter(3) - 8) ! * min(err0 / err, 1.0_r64)

      write(uerr,fmiter) nitert+1, err, 100*ramp

      if (ieee_is_nan(err)) exit relx_corona

      if ((iter > 1 .and. err > err0 * 1.05) .or. (err > 1e4)) then
        write (uerr, '(''diverged: '', Es9.2, '' -> '', Es9.2)') err0, err
        converged = .false.
        exit relx_corona
      end if
      Y(:) = Y + dY * ramp

      where (y_trad < teff * 0.80) y_trad = teff * 0.80
      where (y_temp < teff * 0.80) y_temp = teff * 0.80
      where (.not.ieee_is_normal(y_rho) .or. y_rho < epsilon(1d0)*rhoc) &
            y_rho = epsilon(1d0)*rhoc

      nitert = nitert + 1
      if (cfg_write_all_iters) call saveiter(nitert)

      if (err < 1e-5 .and. err0 / err > 5) then
        write (uerr, '("convergence reached with error = ",ES9.2)') err
        exit relx_corona
      end if

      err0 = err

    end do relx_corona

  end if

  write (uerr,*) '--- DONE'

  !----------------------------------------------------------------------------!
  ! write model data

  write_results: block

    open(33, file = trim(outfn) // '.dat', action = 'write')

    call fillcols(yv, c_, yy)

    do i = 1,ngrid
      write (33,'(I6,*(ES14.5E3))') i, x(i), x(i) / zscale, yy(:,i)
    end do

    close(33)

  end block write_results

  !----------------------------------------------------------------------------!
  ! write some global information

  call wpar_gl(upar)

  write (upar, fmparfc) "alpha", alpha, "alpha parameter"
  write (upar, fmparfc) "eta", zeta, "field rise parameter"
  write (upar, fmparfc) "zeta", zeta, "field rise parameter"
  write (upar, fmparfc) "nu", nu, "reconnection parameter"
  write (upar, fmpare) "radius", radius
  write (upar, fmpare) "zscale", zscale
  write (upar, fmpare) "facc", facc
  write (upar, fmparec) "rho_0", rhoc, "Central density"
  write (upar, fmparec) "temp_0", Tc, "Central gas temperature"
  write (upar, fmpare) "teff", teff
  write (upar, fmpare) "teff_keV", teff * keV_in_kelvin
  write (upar, fmparec) "zdisk_ss73", Hdisk, "Disk height [cm]"
  write (upar, fmparfc) "hdisk_ss73", Hdisk / zscale, "Disk height [cm]"
  write (upar, fmhdr)  "Model information"
  write (upar, fmpari) "model", model
  write (upar, fmpari) "niter", nitert
  write (upar, fmpari) "ngrid", ngrid
  write (upar, fmparl) "converged", converged
  write (upar, fmparl) "has_corona", &
        & (cfg_temperature_method .ne. EQUATION_DIFFUSION)
  write (upar, fmparl) "has_magnetic", .TRUE.
  write (upar, fmparl) "has_conduction", use_conduction

  !----------------------------------------------------------------------------!
  ! save the column information for col2python

  write_columns: block

    open(35, file = trim(outfn) // ".col", action = 'write')

    write(35, fmcol) 'i', 'i4'
    write(35, fmcol) 'z', 'f4'
    write(35, fmcol) 'h', 'f4'

    do i = 1,ncols
      write (35,fmcol) labels(i), 'f4'
    end do

    close(35)

  end block write_columns

  !----------------------------------------------------------------------------!
  ! determine and save the photosphere location

  write_disk_globals: block

    use slf_interpol
    use slf_integrate
    real(dp) :: zphot, ztherm, zhard, zeqbc
    real(dp) :: tavg_therm, tavg_phot, tavg_hard, tavg_eqbc
    real(dp) :: taues_eqbc, taues_therm
    real(dp) :: compfr_phot, compfr_therm
    real(dp) :: compy_phot, compy_therm, compy_hard, compy_eqbc
    real(dp) :: coldens, diskscale
    real(dp), dimension(ngrid) :: compfr

    write (upar,fmhdr) 'Photosphere'

    call interpol(yy(c_tau,:), x, 1d0, zphot)
    write (upar,fmparec) 'zphot', zphot, 'altitude of the photosphere (tau=1)'
    write (upar,fmparf)  'hphot', zphot / zscale

    call interpol(yy(c_tau,:), x, 0.5d0, zhard)
    write (upar,fmparec) 'zhard', zhard, 'altitude of the hard corona (tau=0.5)'
    write (upar,fmparf)  'hhard', zhard / zscale

    call interpol(yy(c_tauth,:), x, 1d0, ztherm)
    write (upar,fmparec) 'ztherm', ztherm, 'thermalization depth (tauth=1)'
    write (upar,fmparf)  'htherm', ztherm / zscale

    call interpol(yy(c_tauth,:), yy(c_taues,:), 1d0, taues_therm)
    write (upar,fmparf)  'taues_therm', taues_therm

    !--------------------------------------------------------------------------!
    ! compton / brehmstrahlung ratio

    compfr = yy(c_coolc,:) / (yy(c_coolc,:) + yy(c_coolb,:))

    call interpol(x, compfr, zphot, compfr_phot)
    write (upar, fmpare) 'compfr_phot', compfr_phot

    call interpol(x, compfr, ztherm, compfr_therm)
    write (upar, fmpare) 'compfr_therm', compfr_therm

    call tabzero(x(ngrid:1:-1), compfr(ngrid:1:-1) - 0.5d0, zeqbc)
    write (upar, fmpare) 'zeqbc', zeqbc
    write (upar, fmpare) 'heqbc', zeqbc / zscale

    call interpol(x, yy(c_taues,:), zeqbc, taues_eqbc)
    write (upar, fmpare) 'taues_eqbc', taues_eqbc

    !--------------------------------------------------------------------------!
    ! average tempearture

    call interpol(x, yy(c_tavg,:), zphot, tavg_phot)
    write (upar, fmparec) 'tavg_phot', tavg_phot, &
          'average temperature up to the photosphere'
    write (upar, fmparf) 'tavg_phot_keV', tavg_phot * keV_in_kelvin

    call interpol(x, yy(c_tavg,:), ztherm, tavg_therm)
    write (upar, fmparec) 'tavg_therm', tavg_therm, &
          'average temperature up to the termalization depth'
    write (upar, fmparf) 'tavg_therm_keV', tavg_therm * keV_in_kelvin

    call interpol(x, yy(c_tavg,:), zhard, tavg_hard)
    write (upar, fmparec) 'tavg_hard', tavg_hard, &
          'average temperature in hot corona'
    write (upar, fmparf) 'tavg_hard_keV', tavg_hard * keV_in_kelvin

    call interpol(x, yy(c_tavg,:), zeqbc, tavg_eqbc)
    write (upar, fmpare) 'tavg_eqbc', tavg_eqbc
    write (upar, fmparf) 'tavg_eqbc_keV', tavg_eqbc * keV_in_kelvin

    !--------------------------------------------------------------------------!
    ! compton Y parameter

    write (upar, fmparfc) 'compy_0', yy(c_compy,1), 'total Y'

    call interpol(x, yy(c_compy,:), zphot, compy_phot)
    write (upar, fmparf) 'compy_phot', compy_phot

    call interpol(x, yy(c_compy,:), ztherm, compy_therm)
    write (upar, fmparf) 'compy_therm', compy_therm

    call interpol(x, yy(c_compy,:), zhard, compy_hard)
    write (upar, fmparf) 'compy_hard', compy_hard

    call interpol(x, yy(c_compy,:), zeqbc, compy_eqbc)
    write (upar, fmparf) 'compy_eqbc', compy_eqbc

    !--------------------------------------------------------------------------!
    ! compute the vertical disk scale and save it

    write (upar, fmhdr)  "column density and disk vertical scale"

    coldens = yy(c_coldens,ngrid)
    write (upar, fmparec) 'coldens', coldens, 'column density'

    diskscale = sqrt(integrate(yy(c_rho,:) * x**2, x) / coldens)

    write (upar, fmparec) 'zdisk', diskscale, &
          'disk vertical scale (second moment)'
    write (upar, fmparfc) 'hdisk', diskscale / zscale, &
          'same but in scale heights'

    !--------------------------------------------------------------------------!
    ! escaping flux - magnetic vs. radiative

    write (upar, fmparf) 'fbfrac_top', &
          yy(c_fmag, ngrid) / (yy(c_frad, ngrid) + yy(c_fmag, ngrid))

    ! some values on the equator
    write (upar, fmpare) 'betarad_0', yy(c_pgas,1) / yy(c_prad,1)

    !--------------------------------------------------------------------------!
    ! evaluate and write the coronal properties

    if (cfg_temperature_method .ne. EQUATION_DIFFUSION .or. cfg_post_corona) then
      write (upar, fmhdr)  "Corona properties"
      block

        real(dp) :: zcor, taucor, frad_disk, tcor, betacor, coldens_disk, fmag_disk, compy_tmin, compfr_tmin

        ! search for temperature minimum
        call findtempmin(x,y_temp,zcor)
        write (upar,fmparec) "zcor", zcor, "height of temperature minimum [cm]"
        write (upar,fmparfc) "hcor", zcor / zscale, "height of temperature minimum [H]"
        call interpol(x, yy(c_temp,:), zcor, tcor)
        write (upar,fmparec) 'tmin', tcor, 'temperature minimum'
        write (upar,fmparec) 'tmin_keV', tcor * keV_in_kelvin, 'same but in keV'

        ! average temperature of the corona
        call interpol(x, yy(c_tavg,:), zcor, tcor)
        write (upar,fmparec) 'tcor', tcor, 'average temperature in corona'
        write (upar,fmparec) 'tcor_keV', tcor * keV_in_kelvin, 'same but in keV'

        ! optical depth of the corona
        call interpol(x, yy(c_tau,:), zcor, taucor)
        write (upar,fmparf) "tau_cor", taucor
        call interpol(x, yy(c_taues,:), zcor, taucor)
        write (upar,fmparf) "taues_cor", taucor
        call interpol(x, yy(c_tauth,:), zcor, taucor)
        write (upar,fmparf) "tauth_cor", taucor

        call interpol(x, yy(c_compy,:), zcor, compy_tmin)
        write (upar,fmparf) 'compy_tmin', compy_tmin

        call interpol(x, compfr, zcor, compfr_tmin)
        write (upar,fmparec) 'compfr_tmin', compfr_tmin, 'compton / brehms. ratio'

        ! energy released in the corona
        call interpol(x,y_frad,zcor,frad_disk)
        write (upar,fmpare) "frad_disk", frad_disk
        write (upar,fmpare) "frad_cor", y_frad(ngrid) - frad_disk
        write (upar,fmparfc) "chicor", 1 - frad_disk / y_frad(ngrid), &
            & "relative amount of radiative flux released in the corona"

        call interpol(x, yy(c_fmag,:), zcor, fmag_disk)
        write (upar, fmparf) 'fbfrac_cor', fmag_disk / (frad_disk + fmag_disk)

        ! magnetic beta
        call interpol(x, yy(c_beta,:), zcor, betacor)
        write (upar,fmparec) 'betacor', betacor, &
            & 'magnetic beta in temperature minimum'

        ! column density: disk vs corona
        call interpol(x, yy(c_coldens,:), zcor, coldens_disk)
        write (upar, fmparec) "coldens_disk", coldens_disk, &
              "column density (disk only)"
        write (upar, fmparec) "coldens_cor", coldens - coldens_disk, &
              "column density (corona only)"
        write (upar, fmparec) "mfrac_cor", (coldens - coldens_disk) / coldens, &
              "mass fraction in the corona"

      end block
    end if
  end block write_disk_globals

  !----------------------------------------------------------------------------!
  ! clean up

  close(upar)

  deallocate(x,x0,Y,M,dY,tau,heat,errmask,ipiv,yy)
  if (allocated(MB)) deallocate(MB)

  !----------------------------------------------------------------------------!

contains

  !----------------------------------------------------------------------------!

  subroutine saveiter(iter)
    use relaxutils
    integer, intent(in) :: iter
    character(256) :: fn
    integer :: i

    write (fn,'(A,".",I0.3,".dat")') trim(outfn),iter
    open(33, file = trim(fn), action = 'write')

    call fillcols(yv, c_, yy)

    writeresults : do i = 1,ngrid
      write (33,'(I6,*(ES14.5E3))') i, x(i), x(i) / zscale, yy(:,i)
    end do writeresults

    close(33)

  end subroutine

  !----------------------------------------------------------------------------!

  subroutine fillcols(yv,c_,yy)
    procedure(funout_t), pointer :: fout
    real(dp) :: kabs,ksct,rhom,tempm,tradm,dx
    real(dp), dimension(:,:), intent(in) :: yv
    integer, dimension(:), intent(in) :: c_
    real(dp), dimension(:,:), intent(inout) :: yy
    integer :: i

    ! select the function appropriate for this model
    call mrx_sel_fout(model, fout)

    do i = 1,ngrid
      call fout(x(i), yv(:,i), yy(:,i))
    end do

    if ( cfg_temperature_method /= EQUATION_BALANCE &
            .and. cfg_post_corona ) then
      block
        use heatbalance, only: heatbil2
        real(dp) :: temp_old
        do i = 1,ngrid
          temp_old = yy(c_temp,i)
          call heatbil2(yy(c_rho,i), yy(c_temp,i), yy(c_trad,i), &
                yy(c_heat,i), .false.)
          yy(c_pgas,i) = yy(c_pgas,i) * yy(c_temp,i) / temp_old
        end do
      end block
    end if

    forall (i = 1:ngrid)
      yy(c_ksct,i) = fksct(yy(c_rho,i), yy(c_temp,i))
      yy(c_kabs,i) = fkabs(yy(c_rho,i), yy(c_temp,i))
      yy(c_kabp,i) = fkabp(yy(c_rho,i), yy(c_temp,i))
      yy(c_kcnd,i) = fkcnd(yy(c_rho,i), yy(c_temp,i))
    end forall

    yy(c_cool0,:) = 4 * cgs_stef * yy(c_rho,:) * yy(c_trad,:)**3 * (yy(c_temp,:) - yy(c_trad,:))
    yy(c_coolb,:) = yy(c_kabp,:) * (1 + (yy(c_temp,:) / yy(c_trad,:))   )  &
                                 * (1 + (yy(c_temp,:) / yy(c_trad,:))**2)
    yy(c_coolc,:) = 4 * yy(c_ksct,:) * cgs_k_over_mec2 * yy(c_trad,:)

    ! compute the optical depths and averaged temperature
    yy(c_tau,ngrid) = 0
    yy(c_taues,ngrid) = 0
    yy(c_tauth,ngrid) = 0
    yy(c_tavg,ngrid) = 0
    yy(c_compy,ngrid) = 0

    integrate_tau : do i = ngrid-1,1,-1
      rhom = (yy(c_rho,i) + yy(c_rho,i+1)) / 2
      if (rhom < 0) rhom = 0
      tempm = (yy(c_temp,i) + yy(c_temp,i+1)) / 2
      tradm = (yy(c_trad,i) + yy(c_trad,i+1)) / 2
      dx = x(i+1) - x(i)

      kabs = merge(fkabs(rhom,tempm), 0.0_dp, tempm > 0)
      ksct = merge(fksct(rhom,tempm), 0.0_dp, tempm > 0)

      yy(c_tau,  i) = yy(c_tau,  i+1) + dx * rhom * (kabs + ksct)
      yy(c_taues,i) = yy(c_taues,i+1) + dx * rhom * ksct
      yy(c_tauth,i) = yy(c_tauth,i+1) + dx * rhom * sqrt(kabs * (kabs + ksct))

      yy(c_tavg, i) = yy(c_tavg, i+1) + dx * rhom * (kabs + ksct) * tempm

      yy(c_compy, i) = yy(c_compy, i+1) + dx * rhom * ksct &
            * 4 * cgs_k_over_mec2 * (tempm - tradm)
    end do integrate_tau

    ! average tempearture
    yy(c_tavg,:) = yy(c_tavg,:) / yy(c_tau,:)
    yy(c_beta,:) = yy(c_pgas,:) / yy(c_pmag,:)
    yy(c_tavg,ngrid) = yy(c_temp,ngrid)

    ! column density
    yy(c_coldens,1) = 0
    integrate_coldens: do i = 1, ngrid-1
      yy(c_coldens, i+1) = yy(c_coldens,i) &
      + (yy(c_rho,i) + yy(c_rho,i+1)) * (x(i+1) - x(i)) / 2
    end do integrate_coldens

  end subroutine

  !----------------------------------------------------------------------------!
  ! searches for temperature minimum in a given array

  subroutine findtempmin(x,temp,xtmin)
    real(dp), intent(in), dimension(:) :: x,temp
    real(dp), intent(out) :: xtmin
    real(dp), dimension(size(x) - 1) :: dtemp,xm
    integer :: i

    do i = 1, size(x)-1
      xm(i) = (x(i) + x(i+1)) / 2
      dtemp(i) = temp(i+1) - temp(i)
    end do

    xtmin = 0
    search_for_minimum : do i = 1,size(dtemp) - 1
      if (dtemp(i) .le. 0 .and. dtemp(i+1) .ge. 0) then
        xtmin = (dtemp(i+1)*xm(i) - dtemp(i)*xm(i+1)) / (dtemp(i+1) - dtemp(i))
        exit search_for_minimum
      end if
    end do search_for_minimum

  end subroutine

  !----------------------------------------------------------------------------!
  ! searches for zero in the array

  subroutine tabzero(x,y,x0)
    real(dp), intent(in), dimension(:) :: x,y
    real(dp), intent(inout) :: x0
    integer :: i

    if (size(x) /= size(y)) error stop

    search_for_zero : do i = 1, size(y)-1
      if (y(i) * y(i+1) .le. 0) then
        x0 = (y(i+1)*x(i) - y(i)*x(i+1)) / (y(i+1) - y(i))
        exit search_for_zero
      end if
    end do search_for_zero

  end subroutine

  !----------------------------------------------------------------------------!

  subroutine rdconf(cfg)
    integer :: errno
    type(config), intent(inout) :: cfg
    character(len=2048) :: buf

    call mincf_get(cfg, "alpha", buf, errno)
    if ( iand(errno, mincf_not_found) .ne. 0 )  then
      error stop "Magnetic alpha-parameter (key: alpha) is REQUIRED!"
    end if
    read (buf,*) alpha

    call mincf_get(cfg, "radius", buf, errno)
    if ( iand(errno, mincf_not_found) .ne. 0 )  then
      error stop "Radius (key: radius) is REQUIRED!"
    end if
    read (buf,*) radius

    call mincf_get(cfg, "zeta", buf, errno)
    if ( iand(errno, mincf_not_found) .ne. 0 )  then
      error stop "Magnetic field parameter (key: zeta) is REQUIRED!"
    end if
    read (buf,*) zeta

    call mincf_get(cfg, "nu", buf, errno)
    if ( iand(errno, mincf_not_found) .eq. 0 )  then
      read (buf,*) nu
    else
      write (0,*) "warning: assuming nu = 0"
      nu = 0
    end if

  end subroutine

  !----------------------------------------------------------------------------!

end program
