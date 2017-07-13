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
  integer, dimension(3) :: niter = [ 36, 8, 52 ]
  real(dp), allocatable, target :: x(:), x0(:), Y(:), dY(:), M(:,:)
  logical, dimension(:), allocatable :: errmask
  real(dp), pointer, dimension(:) :: y_rho, y_temp, y_frad, y_pmag, y_Trad
  real(dp), allocatable, dimension(:) :: tau, d_frad
  real(dp) :: rhoc, Tc, Hdisk, err, err0
  logical :: user_ff, user_bf
  integer, dimension(6) :: c_
  integer, parameter :: upar = 92

  !----------------------------------------------------------------------------!

  ngrid = 720
  htop = 48

  !----------------------------------------------------------------------------!

  call rdargvgl
  call rdargvrx
  call mincf_read(cfg)
  call rdconfgl(cfg)
  call rdconf(cfg)
  call mincf_free(cfg)

  open(upar, file = trim(outfn) // '.txt', action = 'write')

  user_ff = use_opacity_ff
  user_bf = use_opacity_bf

  !----------------------------------------------------------------------------!
  ! calculate the global parameters
  call mrx_init(mbh, mdot, radius, alpha)

  ! get the model number
  model = mrx_number( 'D', .TRUE., .FALSE. )
  ny = mrx_ny(model)
  call mrx_sel_hash(model, C_)

  allocate( x(ngrid), x0(ngrid), Y(ny*ngrid), dY(ny*ngrid),  M(ny*ngrid,ny*ngrid), tau(ngrid), d_frad(ngrid) )
  allocate(errmask(ny*ngrid))

  !----------------------------------------------------------------------------!

  call apx_estim(mbh, mdot, radius, alpha, rhoc, Tc, Hdisk)
  call apx_refin(mbh, mdot, radius, alpha, rhoc, Tc, Hdisk)

  !----------------------------------------------------------------------------!

  forall (i = 1:ngrid)
    x(i) = space_linear(i,ngrid,htop) * Hdisk
    x0(i) = space_linear(i,ngrid,htop) / htop
  end forall

  y_rho   => Y(C_(1)::ny)
  y_temp  => Y(C_(2)::ny)
  y_trad  => Y(C_(3)::ny)
  y_frad  => Y(C_(4)::ny)
  y_pmag  => Y(C_(5)::ny)

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

    call mrx_advance(model, x, Y, M, dY, errno)

    forall (i = 1:ngrid*ny) errmask(i) = (Y(i) .ne. 0) &
        & .and. ieee_is_normal(dY(i))
    err = sqrt(sum((dY/Y)**2, errmask) / count(errmask))
    write(uerr,'(I5,ES9.2)') nitert+1, err

    Y(:) = Y + dY * ramp3(iter,niter(1))

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

      call mrx_advance(model, x, Y, M, dY, errno)

      forall (i = 1:ngrid*ny) errmask(i) = (Y(i) .ne. 0) &
          & .and. ieee_is_normal(dY(i))
      err = sqrt(sum((dY/Y)**2, errmask) / count(errmask))
      write(uerr,'(I5,ES9.2)') nitert+1, err

      Y(:) = Y + dY * ramp2(iter,niter(2),0.5d0)

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

  if ( cfg_temperature_method /= EQUATION_DIFFUSION ) then

    write (uerr,*) '--- corona is on'

    err0 = 0

    call mrx_transfer(model, mrx_number(cfg_temperature_method, .TRUE., .FALSE.), x, Y)

    call mrx_sel_ny(model,ny)
    call mrx_sel_hash(model,c_)

    deallocate(dY,M,errmask)
    allocate(dY(ny*ngrid), M(ny*ngrid,ny*ngrid),errmask(ny*ngrid))

    y_rho   => Y(C_(1)::ny)
    y_temp  => Y(C_(2)::ny)
    y_trad  => Y(C_(3)::ny)
    y_frad  => Y(C_(4)::ny)
    y_pmag  => Y(C_(5)::ny)

    call deriv(x, y_frad, d_frad)

    forall (i = 1:ngrid)
      y_temp(i) = d_frad(i) * (cgs_mel * cgs_c**2) &
            & / (16 * cgs_boltz * (cgs_kapes*y_rho(i)) &
            & * (cgs_stef*y_trad(i)**4) )
      y_temp(i) = sqrt(y_temp(i)**2 + y_trad(i)**2)
    end forall

    nitert = nitert + 1
    if (cfg_write_all_iters) call saveiter(nitert)

    relx_corona : do iter = 1,niter(3)

      call mrx_advance(model, x, Y, M, dY, errno)

      if (errno .ne. 0) exit relx_corona

      forall (i = 1:ngrid*ny) errmask(i) = (Y(i) .ne. 0) &
          & .and. ieee_is_normal(dY(i))
      err = sqrt(sum((dY/Y)**2, errmask) / count(errmask))
      write(uerr,'(I5,ES9.2)') nitert+1, err

      if (ieee_is_nan(err)) exit relx_corona

      Y(:) = Y + dY * ramp3(iter, niter(3))

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

  open(33, file = trim(outfn) // '.dat', action = 'write')
  call saveresult(33)
  close(33)

  !----------------------------------------------------------------------------!

  call wpar_gl(upar)
  write (upar, fmparfc) "alpha", alpha, "Alpha parameter"
  write (upar, fmparfc) "zeta", zeta, "Zeta parameter"
  write (upar, fmpare) "radius", radius
  write (upar, fmpare) "zscale", fzscale(mbh,mdot,radius)
  write (upar, fmparec) "rho_0", rhoc, "Central density"
  write (upar, fmparec) "temp_0", Tc, "Central gas temperature"
  write (upar, fmparec) "Hdisk", Hdisk, "Disk height [cm]"
  write (upar, fmhdr)  "Model information"
  write (upar, fmpari) "model", model
  write (upar, fmpari) "niter", nitert
  write (upar, fmparl) "has_corona", &
        & (cfg_temperature_method == EQUATION_BALANCE)
  write (upar, fmparl) "has_magnetic", .TRUE.
  write (upar, fmparl) "has_conduction", .FALSE.

  open(35, file = trim(outfn) // ".col", action = 'write')
  write(35, fmcol) 'i', 'i4'
  write(35, fmcol) 'z', 'f4'
  write(35, fmcol) 'h', 'f4'
  write(35, fmcol) 'tau', 'f4'
  write(35, fmcol) 'rho', 'f4'
  write(35, fmcol) 'temp', 'f4'
  write(35, fmcol) 'trad', 'f4'
  write(35, fmcol) 'frad', 'f4'
  write(35, fmcol) 'ksct', 'f4'
  write(35, fmcol) 'kabs', 'f4'
  write(35, fmcol) 'kcnd', 'f4'
  write(35, fmcol) 'pgas', 'f4'
  write(35, fmcol) 'prad', 'f4'
  write(35, fmcol) 'pmag', 'f4'
  write(35, fmcol) 'beta', 'f4'
  write(35, fmcol) 'd*frad', 'f4'
  close(35)

  close(upar)

  deallocate(x,x0,Y,M,dY,tau,d_frad,errmask)

contains

  !----------------------------------------------------------------------------!
  subroutine saveiter(iter)
    use relaxutils
    integer, intent(in) :: iter
    character(256) :: fn
    write (fn,'(A,".",I0.3,".dat")') trim(outfn),iter
    open(33, file = trim(fn), action = 'write')
    call saveresult(33)
    close(33)
  end subroutine

  subroutine saveresult(u)
    integer, intent(in) :: u
    real(r64) :: pgas, prad
    integer :: i
    call mktaues(x,y_rho,y_temp,tau)
    call deriv(x, y_frad, d_frad)

    do i = 1,ngrid
      pgas =  cgs_k_over_mh / miu * y_rho(i) * y_temp(i)
      prad = cgs_a / 3 * y_trad(i)**4
      write (u,'(I6,15ES14.5E3)') i, x(i), x(i) / zscale, tau(i), &
      y_rho(i), y_temp(i), y_trad(i), y_frad(i),  &
      fksct(y_rho(i), y_temp(i)), fkabs(y_rho(i), y_temp(i)), &
      fkcnd(y_rho(i), y_temp(i)), &
      pgas, prad, y_pmag(i), pgas / y_pmag(i), d_frad(i)
    end do
  end subroutine

  !----------------------------------------------------------------------------!
  subroutine rdconf(cfg)
    integer :: errno
    type(config), intent(inout) :: cfg
    character(len=2048) :: buf

    call mincf_get(cfg, "alpha", buf)
    if ( mincf_failed(cfg) )  then
      error stop "Magnetic alpha-parameter (key: alpha) is REQUIRED!"
    end if
    read (buf,*) alpha

    call mincf_get(cfg, "radius", buf)
    if ( mincf_failed(cfg) )  then
      error stop "Radius (key: radius) is REQUIRED!"
    end if
    read (buf,*) radius

    call mincf_get(cfg, "zeta", buf)
    if ( mincf_failed(cfg) )  then
      error stop "Magnetic field parameter (key: zeta) is REQUIRED!"
    end if
    read (buf,*) zeta

  end subroutine

end program
