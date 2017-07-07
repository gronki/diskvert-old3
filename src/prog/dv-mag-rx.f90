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
  use ss73solution, only: apxdisk, apx_estim, apx_refine
  use grid, only: space_linear, space_linlog

  !----------------------------------------------------------------------------!
  implicit none

  type(config) :: cfg
  integer :: model, errno
  integer :: ny = 3, i, iter, nitert = 0
  integer, dimension(3) :: niter = [ 36, 8, 256 ]
  character(2**8) :: fn
  real(dp), allocatable, target :: x(:), x0(:), Y(:), dY(:), M(:,:)
  real(dp), pointer, dimension(:) :: y_rho, y_temp, y_frad, y_pmag, y_Trad
  real(dp), allocatable, dimension(:) :: tau, d_frad
  real(dp) :: rhoc, Tc, Hdisk, err, t
  logical :: user_ff, user_bf
  integer, dimension(6) :: c_

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

  user_ff = use_opacity_ff
  user_bf = use_opacity_bf

  !----------------------------------------------------------------------------!
  ! calculate the global parameters
  call mrx_init(mbh, mdot, radius, alpha)

  ! get the model number
  model = mrx_number( .FALSE., .TRUE., .FALSE. )
  ny = mrx_ny(model)
  call mrx_sel_hash(model, C_)

  allocate( x(ngrid), x0(ngrid), Y(ny*ngrid), dY(ny*ngrid), M(ny*ngrid,ny*ngrid), tau(ngrid), d_frad(ngrid) )

  !----------------------------------------------------------------------------!

  call apx_estim(mbh, mdot, radius, alpha, rhoc, Tc, Hdisk)
  call apx_refine(mbh, mdot, radius, alpha, rhoc, Tc, Hdisk)

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

  y_frad = x0 * facc
  y_temp =exp(-0.5*(x/Hdisk)**2) * (Tc - 0.841 * Teff) + 0.841 * Teff
  y_rho =  rhoc * (exp(-0.5*(x/Hdisk)**2) + epsilon(1.0_dp))

  associate (beta_0 => 2 * zeta / alpha - 1)
    y_pmag = 2 * cgs_k_over_mh * y_rho * y_temp   &
          & / (beta_0 * exp(- 0.5 * (x / Hdisk)**2 ) + 2e-2)
  end associate

  !----------------------------------------------------------------------------!

  if (cfg_write_all_iters) call saveiter(0)

  use_opacity_ff = .false.
  use_opacity_bf = .false.
  relx_opacity_es : do iter = 1,niter(1)

    call mrx_advance(model, x, Y, M, dY, errno)

    err = sum(merge((dY/Y)**2, 0d0, Y.ne.0)) / sum(merge(1, 0, Y.ne.0))
    write(ulog,'(I5,Es12.4)') nitert+1, err

    Y = Y + dY * ramp3(iter,niter(1))
    ! y_rho = merge(y_rho, epsilon(1d0)*rhoc, y_rho > epsilon(1d0)*rhoc)
    ! y_temp = merge(y_temp, teff / 2, y_temp > teff / 2)

    nitert = nitert + 1
    if (cfg_write_all_iters) call saveiter(nitert)

  end do relx_opacity_es


  !----------------------------------------------------------------------------!

  use_opacity_ff = user_ff
  use_opacity_bf = user_bf

  if ( use_opacity_bf .or. use_opacity_ff ) then

    write (ulog,*) '--- ff+bf opacity is on'

    relx_opacity_full : do iter = 1,niter(2)

      call mrx_advance(model, x, Y, M, dY, errno)

      err = sum(merge((dY/Y)**2, 0d0, Y.ne.0)) / sum(merge(1, 0, Y.ne.0))
      write(ulog,'(I5,Es12.4)') nitert+1, err

      Y = Y + dY * ramp2(iter,niter(2),0.5d0)
      ! y_rho = merge(y_rho, epsilon(1d0)*rhoc, y_rho > epsilon(1d0)*rhoc)
      ! y_temp = merge(y_temp, teff / 2, y_temp > teff / 2)

      nitert = nitert + 1
      if (cfg_write_all_iters) call saveiter(nitert)

    end do relx_opacity_full


  end if

  !----------------------------------------------------------------------------!

  if ( cfg_temperature_method == EQUATION_BALANCE ) then

    write (ulog,*) '--- corona is on'

    call mrx_transfer(model, mrx_number( .TRUE., .TRUE., .FALSE. ), x, Y)

    ny = mrx_ny(model)
    call mrx_sel_hash(model,c_)

    deallocate(dY,M)
    allocate(dY(ny*ngrid), M(ny*ngrid,ny*ngrid))

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

      err = sum(merge((dY/Y)**2, 0d0, Y.ne.0)) / sum(merge(1, 0, Y.ne.0))
      write(ulog,'(I5,Es12.4)') nitert+1, err

      if (ieee_is_nan(err)) exit relx_corona

      Y = Y + dY * ramp3(iter, niter(3))

      where (y_trad < teff * 0.80) y_trad = teff * 0.80
      where (y_temp < teff * 0.80) y_temp = teff * 0.80
      where (.not.ieee_is_normal(y_rho) .or. y_rho < epsilon(1d0)*rhoc) &
            y_rho = epsilon(1d0)*rhoc

      nitert = nitert + 1
      if (cfg_write_all_iters) call saveiter(nitert)

    end do relx_corona

  end if

  write (ulog,*) '--- DONE'

  open(33, file = trim(outfn) // '.dat', action = 'write')
  call saveresult(33)
  close(33)

  !----------------------------------------------------------------------------!

  open(newunit = upar, file = trim(outfn) // '.txt', action = 'write')
  call wpar_gl(upar)
  write (upar, fmparfc) "alpha", alpha, "Alpha parameter"
  write (upar, fmparfc) "zeta", zeta, "Zeta parameter"
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
  close(upar)

  open(34, file = trim(outfn) // ".col", action = 'write')
  write(34, fmcol) 'i', 'i4'
  write(34, fmcol) 'z', 'f4'
  write(34, fmcol) 'h', 'f4'
  write(34, fmcol) 'tau', 'f4'
  write(34, fmcol) 'rho', 'f4'
  write(34, fmcol) 'temp', 'f4'
  write(34, fmcol) 'trad', 'f4'
  write(34, fmcol) 'frad', 'f4'
  write(34, fmcol) 'ksct', 'f4'
  write(34, fmcol) 'kabs', 'f4'
  write(34, fmcol) 'kcnd', 'f4'
  write(34, fmcol) 'pgas', 'f4'
  write(34, fmcol) 'prad', 'f4'
  write(34, fmcol) 'pmag', 'f4'
  write(34, fmcol) 'beta', 'f4'
  write(34, fmcol) 'd*frad', 'f4'
  close(34)

  deallocate(x,x0,Y,M,dY,tau,d_frad)

contains

  !----------------------------------------------------------------------------!
  subroutine saveiter(iter)
    use relaxutils
    integer, intent(in) :: iter
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
