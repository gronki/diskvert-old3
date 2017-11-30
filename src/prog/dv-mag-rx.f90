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
  use grid, only: space_linear, space_linlog

  !----------------------------------------------------------------------------!
  implicit none

  type(config) :: cfg
  integer :: model, errno
  integer :: ny = 3, i, iter, nitert = 0
  integer, dimension(3) :: niter = [ 24, 8, 72 ]
  real(dp), allocatable, target :: x(:), x0(:), Y(:), dY(:), M(:,:), YY(:,:)
  real(dp), pointer :: yv(:,:)
  integer, dimension(:), allocatable :: ipiv
  logical, dimension(:), allocatable :: errmask
  real(dp), pointer, dimension(:) :: y_rho, y_temp, y_frad, y_pmag, y_Trad
  real(dp), allocatable, dimension(:) :: tau, heat
  real(dp) :: rhoc, Tc, Hdisk, err, err0, ramp
  character(*), parameter :: fmiter = '(I5,2X,ES9.2,2X,F5.1,"%")'
  logical :: user_ff, user_bf
  integer, dimension(6) :: c_
  integer, parameter :: upar = 92

  integer, parameter :: ncols  = 19, &
      c_rho = 1, c_temp = 2, c_trad = 3, &
      c_pgas = 4, c_prad = 5, c_pmag = 6, &
      c_frad = 7, c_fmag = 8, c_fcnd = 9, &
      c_heat = 10, c_vrise = 11, &
      c_ksct = 12, c_kabs = 13, c_kcnd = 14, &
      c_tau = 15, c_taues = 16, c_tauth = 17, &
      c_tavg = 18, c_beta = 19
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
  labels(c_kcnd) = 'kcnd'
  labels(c_tau) = 'tau'
  labels(c_taues) = 'taues'
  labels(c_tauth) = 'tauth'
  labels(c_tavg) = 'tavg'
  labels(c_beta) = 'beta'

  !----------------------------------------------------------------------------!

  ngrid = -1
  htop = 120

  !----------------------------------------------------------------------------!

  call rdargvgl
  call rdargvrx
  call mincf_read(cfg)
  call rdconfgl(cfg)
  call rdconf(cfg)
  call mincf_free(cfg)

  if (ngrid .eq. -1) ngrid = ceiling(8*htop)

  open(upar, file = trim(outfn) // '.txt', action = 'write')

  user_ff = use_opacity_ff
  user_bf = use_opacity_bf

  !----------------------------------------------------------------------------!
  ! calculate the global parameters
  call mrx_init(mbh, mdot, radius, alpha)

  ! get the model number
  model = mrx_number( 'D', .TRUE., .FALSE. )
  call mrx_sel_nvar(model, ny)
  call mrx_sel_hash(model, C_)

  allocate( x(ngrid), x0(ngrid), Y(ny*ngrid), dY(ny*ngrid),  M(ny*ngrid,ny*ngrid), tau(ngrid), heat(ngrid), YY(ncols,ngrid) )
  allocate(errmask(ny*ngrid), ipiv(ny*ngrid))

  !----------------------------------------------------------------------------!

  call apx_estim(mbh, mdot, radius, alpha, rhoc, Tc, Hdisk)
  call apx_refin(mbh, mdot, radius, alpha, rhoc, Tc, Hdisk)

  !----------------------------------------------------------------------------!

  forall (i = 1:ngrid)
    x(i) = space_linear(i,ngrid,htop) * zscale
    x0(i) = space_linear(i,ngrid,htop) / htop
  end forall

  y_rho   => Y(C_(1)::ny)
  y_temp  => Y(C_(2)::ny)
  y_trad  => Y(C_(3)::ny)
  y_frad  => Y(C_(4)::ny)
  y_pmag  => Y(C_(5)::ny)
  YV(1:ny,1:ngrid) => Y

  !----------------------------------------------------------------------------!

  associate (beta_0 => 2 * zeta / alpha - 1)
    if (beta_0 < 0) error stop "MAGNETIC BETA cannot be negative! " &
          & // "zeta > alpha / 2!!!"
    write (upar,fmpare) 'beta_0', beta_0
    forall (i = 1:ngrid)
      y_frad(i) = x0(i) * facc
      y_temp(i) = (1 - x0(i)) * (Tc - 0.841 * Teff) + 0.841 * Teff
      y_rho(i) =  rhoc * (exp(-0.5*(x(i)/Hdisk)**2) + 1e-8)

      y_pmag(i) = 2 * cgs_k_over_mh * y_rho(i) * y_temp(i)   &
            & / (beta_0 * exp(- 0.5 * (x(i) / Hdisk)**2 ) + 1e-2)
    end forall
  end associate

  !----------------------------------------------------------------------------!

  if (cfg_write_all_iters) call saveiter(0)

  use_opacity_ff = .false.
  use_opacity_bf = .false.
  err0 = 0

  relx_opacity_es : do iter = 1,niter(1)

    call mrx_matrix(model, x, Y, M, dY)
    call dgesv(size(M,2), 1, M, size(M,1), ipiv, dY, size(dY), errno)

    forall (i = 1:ngrid*ny) errmask(i) = (Y(i) .ne. 0) &
        & .and. ieee_is_normal(dY(i))
    err = sqrt(sum((dY/Y)**2, errmask) / count(errmask))
    ramp = ramp4(iter, niter(1), 5d-2, 0.25d0)
    write(uerr,fmiter) nitert+1, err, 100*ramp

    Y(:) = Y + dY * ramp

    nitert = nitert + 1
    if (cfg_write_all_iters) call saveiter(nitert)

    if (err < 1e-6 .and. err0 / err > 10) then
      write (uerr, '("convergence reached with error = ",ES9.2)') err
      exit relx_opacity_es
    end if

    err0 = err

  end do relx_opacity_es


  !----------------------------------------------------------------------------!

  use_opacity_ff = user_ff
  use_opacity_bf = user_bf

  err0 = 0

  if ( use_opacity_bf .or. use_opacity_ff ) then

    write (uerr,*) '--- ff+bf opacity is on'

    relx_opacity_full : do iter = 1,niter(2)

      call mrx_matrix(model, x, Y, M, dY)
      call dgesv(size(M,2), 1, M, size(M,1), ipiv, dY, size(dY), errno)

      forall (i = 1:ngrid*ny) errmask(i) = (Y(i) .ne. 0) &
          & .and. ieee_is_normal(dY(i))
      err = sqrt(sum((dY/Y)**2, errmask) / count(errmask))
      ramp = ramp4(iter, niter(2), 2d-1, 0.5d0)
      write(uerr,fmiter) nitert+1, err, 100*ramp

      Y(:) = Y + dY * ramp

      nitert = nitert + 1
      if (cfg_write_all_iters) call saveiter(nitert)

      if (err < 1e-6 .and. err0 / err > 10) then
        write (uerr, '("convergence reached with error = ",ES9.2)') err
        exit relx_opacity_full
      end if

      err0 = err

    end do relx_opacity_full


  end if

  !----------------------------------------------------------------------------!

  if ( cfg_temperature_method .ne. EQUATION_DIFFUSION ) then

    write (uerr,*) '--- corona is on'

    err0 = 0

    call mrx_transfer(model, &
        mrx_number(cfg_temperature_method, .TRUE., .FALSE.), x, Y)

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

    call deriv(x, y_frad, heat)

    forall (i = 1:ngrid)
      y_temp(i) = heat(i) * (cgs_mel * cgs_c**2) &
            & / (16 * cgs_boltz * (cgs_kapes*y_rho(i)) &
            & * (cgs_stef*y_trad(i)**4) )
      y_temp(i) = sqrt(y_temp(i)**2 + y_trad(i)**2)
    end forall

    nitert = nitert + 1
    if (cfg_write_all_iters) call saveiter(nitert)

    !--------------------------------------------------------------------------!

    relx_corona : do iter = 1,niter(3)

      call mrx_matrix(model, x, Y, M, dY)
      call dgesv(size(M,2), 1, M, size(M,1), ipiv, dY, size(dY), errno)

      if (errno .ne. 0) exit relx_corona

      forall (i = 1:ngrid*ny) errmask(i) = (Y(i) .ne. 0) &
          & .and. ieee_is_normal(dY(i))
      err = sqrt(sum((dY/Y)**2, errmask) / count(errmask))
      ramp = ramp4(iter, niter(3), 1d-2, 0d1)

      write(uerr,fmiter) nitert+1, err, 100*ramp

      if (ieee_is_nan(err)) exit relx_corona

      Y(:) = Y + dY * ramp

      where (y_trad < teff * 0.80) y_trad = teff * 0.80
      where (y_temp < teff * 0.80) y_temp = teff * 0.80
      where (.not.ieee_is_normal(y_rho) .or. y_rho < epsilon(1d0)*rhoc) &
            y_rho = epsilon(1d0)*rhoc

      nitert = nitert + 1
      if (cfg_write_all_iters) call saveiter(nitert)

      if (err < 1e-6 .and. err0 / err > 10) then
        write (uerr, '("convergence reached with error = ",ES9.2)') err
        exit relx_corona
      end if

      err0 = err

    end do relx_corona

  end if

  write (uerr,*) '--- DONE'

  !----------------------------------------------------------------------------!

  open(33, file = trim(outfn) // '.dat', action = 'write')
  call fillcols(yv, c_, yy)
  writeresults : do i = 1,ngrid
    write (33,'(I6,*(ES14.5E3))') i, x(i), x(i) / zscale, yy(:,i)
  end do writeresults
  close(33)

  !----------------------------------------------------------------------------!

  call wpar_gl(upar)
  write (upar, fmparfc) "alpha", alpha, "Alpha parameter"
  write (upar, fmparfc) "zeta", zeta, "Zeta parameter"
  write (upar, fmpare) "radius", radius
  write (upar, fmpare) "zscale", zscale
  write (upar, fmparec) "rho_0", rhoc, "Central density"
  write (upar, fmparec) "temp_0", Tc, "Central gas temperature"
  write (upar, fmpare) "teff", teff
  write (upar, fmparec) "Hdisk", Hdisk, "Disk height [cm]"
  write (upar, fmhdr)  "Model information"
  write (upar, fmpari) "model", model
  write (upar, fmpari) "niter", nitert
  write (upar, fmparl) "has_corona", &
        & (cfg_temperature_method .ne. EQUATION_DIFFUSION)
  write (upar, fmparl) "has_magnetic", .TRUE.
  write (upar, fmparl) "has_conduction", .FALSE.

  !----------------------------------------------------------------------------!

  open(35, file = trim(outfn) // ".col", action = 'write')
  write(35, fmcol) 'i', 'i4'
  write(35, fmcol) 'z', 'f4'
  write(35, fmcol) 'h', 'f4'
  do i = 1,ncols
    write (35,fmcol) labels(i), 'f4'
  end do
  close(35)

  !----------------------------------------------------------------------------!

  if (cfg_temperature_method .ne. EQUATION_DIFFUSION) then
    write (upar, fmhdr)  "Corona properties"
    coronal_properties : block

      use slf_interpol
      real(dp) :: zcor, taucor, frad_disk, tcor, betacor

      ! search for temperature minimum
      call findtempmin(x,y_temp,zcor)
      write (upar,fmparec) "zcor", zcor, "height of temperature minimum [cm]"
      write (upar,fmparfc) "hcor", zcor / zscale, "height of temperature minimum [H]"
      call interpol(x, yy(c_temp,:), zcor, tcor)
      write (upar,fmparec) 'tmin', tcor, 'temperature minimum'

      ! average temperature of the corona
      call interpol(x, yy(c_tavg,:), zcor, tcor)
      write (upar,fmparec) 'tcor', tcor, 'average temperature in corona'

      ! optical depth of the corona
      call interpol(x, yy(c_tau,:), zcor, taucor)
      write (upar,fmparf) "tau_cor", taucor
      call interpol(x, yy(c_taues,:), zcor, taucor)
      write (upar,fmparf) "taues_cor", taucor
      call interpol(x, yy(c_tauth,:), zcor, taucor)
      write (upar,fmparf) "tauth_cor", taucor

      ! energy released in the corona
      call interpol(x,y_frad,zcor,frad_disk)
      write (upar,fmpare) "frad_disk", frad_disk
      write (upar,fmpare) "frad_cor", y_frad(ngrid) - frad_disk
      write (upar,fmparfc) "chicor", 1 - frad_disk / y_frad(ngrid), &
      & "relative amount of radiative flux released in the corona"

      ! magnetic beta
      call interpol(x, yy(c_beta,:), zcor, betacor)
      write (upar,fmparec) 'betacor', betacor, &
        & 'magnetic beta in temperature minimum'

    end block coronal_properties
  end if

  !----------------------------------------------------------------------------!

  PhotosphereLocation : block
    use slf_interpol
    real(dp) :: zphot, ztherm

    write (upar,fmhdr) 'Photosphere'

    call interpol(yy(c_tau,:), x, 1d0, zphot)
    write (upar,fmparec) 'zphot', zphot, 'altitude of the photosphere (tau=1)'

    call interpol(yy(c_tauth,:), x, 1d0, ztherm)
    write (upar,fmparec) 'ztherm', ztherm, 'thermalization depth (tauth=1)'
  end block PhotosphereLocation


  !----------------------------------------------------------------------------!

  close(upar)

  deallocate(x,x0,Y,M,dY,tau,heat,errmask,ipiv,yy)

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
    real(dp) :: kabs,ksct,rhom,tempm,dx
    real(dp), dimension(:,:), intent(in) :: yv
    integer, dimension(:), intent(in) :: c_
    real(dp), dimension(:,:), intent(inout) :: yy
    integer :: i

    ! select the function appropriate for this model
    call mrx_sel_fout(model, fout)

    do i = 1,ngrid
      call fout(x(i), yv(:,i), yy(:,i))
      yy(c_ksct,i) = fksct(yy(c_rho,i),yy(c_temp,i))
      yy(c_kabs,i) = fkabs(yy(c_rho,i),yy(c_temp,i))
      yy(c_kcnd,i) = fkcnd(yy(c_rho,i),yy(c_temp,i))
    end do

    ! compute the optical depths and averaged temperature
    yy(c_tau,ngrid) = 0
    yy(c_taues,ngrid) = 0
    yy(c_tauth,ngrid) = 0
    yy(c_tavg,ngrid) = 0

    integrate_tau : do i = ngrid-1,1,-1
      rhom = (yy(c_rho,i) + yy(c_rho,i+1)) / 2
      if (rhom < 0) rhom = 0
      tempm = (yy(c_temp,i) + yy(c_temp,i+1)) / 2
      dx = x(i+1) - x(i)
      kabs = merge(fkabs(rhom,tempm), 0.0_dp, tempm > 0)
      ksct = merge(fksct(rhom,tempm), 0.0_dp, tempm > 0)
      yy(c_tau,  i) = yy(c_tau,  i+1) + dx * rhom * (kabs + ksct)
      yy(c_tavg, i) = yy(c_tavg, i+1) + dx * rhom * (kabs + ksct) * tempm
      yy(c_taues,i) = yy(c_taues,i+1) + dx * rhom * ksct
      yy(c_tauth,i) = yy(c_tauth,i+1) + dx * rhom * sqrt(kabs*(kabs+ksct))
    end do integrate_tau

    ! average tempearture
    forall (i = 1:ngrid)
      yy(c_tavg,i) = yy(c_tavg,i) / yy(c_tau,i)
      yy(c_beta,i) = yy(c_pgas,i) / yy(c_pmag,i)
    end forall
    yy(c_tavg,ngrid) = yy(c_temp,ngrid)

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

  end subroutine

  !----------------------------------------------------------------------------!

end program
