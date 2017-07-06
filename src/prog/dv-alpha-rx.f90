program dv_alpha_relax

  use confort
  use globals
  use settings
  use iso_fortran_env, only: sp => real32, dp => real64
  use ieee_arithmetic, only: ieee_is_normal, ieee_is_nan
  use relaxutils
  use fileunits
  use relaxation
  use rxsettings
  use ss73solution, only: apxdisk, apx_estim, apx_refine
  use grid, only: space_linear, space_linlog

  !----------------------------------------------------------------------------!
  implicit none

  type(config) :: cfg
  integer :: model
  integer :: ny = 3, i, iter, nitert = 0
  integer, dimension(2) :: niter = [ 36, 12 ]
  character(2**8) :: fn
  real(dp), allocatable, target :: x(:), x0(:), Y(:), dY(:), M(:,:)
  real(dp), pointer, dimension(:) :: y_rho, y_temp, y_frad, tau, y_trad
  real(dp) :: rhoc, Tc, Hdisk, err, t
  logical :: user_ff, user_bf
  integer, dimension(6) :: C_

  !----------------------------------------------------------------------------!

  ngrid = 2**7
  htop = 4

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
  model = mrx_number( .FALSE., .FALSE., .FALSE. )
  ny = mrx_ny(model)
  call mrx_sel_hash(model, C_)

  allocate( x(ngrid), x0(ngrid), Y(ny*ngrid), dY(ny*ngrid), M(ny*ngrid,ny*ngrid), tau(ngrid) )

  !----------------------------------------------------------------------------!

  call apx_estim(mbh, mdot, radius, alpha, rhoc, Tc, Hdisk)
  call apx_refine(mbh, mdot, radius, alpha, rhoc, Tc, Hdisk)

  !----------------------------------------------------------------------------!

  forall (i = 1:ngrid)
    x(i) = space_linlog(i,ngrid,htop) * Hdisk
    x0(i) = space_linlog(i,ngrid,htop) / htop
  end forall

  y_rho   => Y(C_(1)::ny)
  y_temp  => Y(C_(2)::ny)
  y_trad  => Y(C_(2)::ny)
  y_frad  => Y(C_(4)::ny)

  !----------------------------------------------------------------------------!

  y_frad = x0 * facc
  y_temp = exp(-0.5*(x/Hdisk)**2) * (Tc - 0.841 * Teff) + 0.841 * Teff
  y_rho =  rhoc * (exp(-0.5*(x/Hdisk)**2) + 1e-5)

  !----------------------------------------------------------------------------!
  if (cfg_write_all_iters) call saveiter(0)

  use_opacity_ff = .false.
  use_opacity_bf = .false.

  relx_opacity_es : do iter = 1,niter(1)

    call mrx_advance(model, x, Y, M, dY)

    err = sum(merge((dY/Y)**2, 0d0, Y.ne.0)) / sum(merge(1, 0, Y.ne.0))
    write(ulog,'(I5,Es12.4)') iter, err

    Y = Y + dY * ramp1(iter,niter(1))

    y_rho = merge(y_rho, 0.0_dp, ieee_is_normal(y_rho) .and. (y_rho > 0))
    y_temp = merge(y_temp, teff / 2, y_temp > teff / 2)

    if (cfg_write_all_iters) call saveiter(iter)

  end do relx_opacity_es

  nitert = nitert + niter(1)

  !----------------------------------------------------------------------------!

  use_opacity_ff = user_ff
  use_opacity_bf = user_bf

  if ( use_opacity_ff .or. use_opacity_bf ) then
    write (ulog,*) '--- bf+ff opacity is ON'

    relx_opacity_full : do iter = 1,niter(2)

      call mrx_advance(model, x, Y, M, dY)

      err = sum(merge((dY/Y)**2, 0d0, Y.ne.0)) / sum(merge(1, 0, Y.ne.0))
      write(ulog,'(I5,Es12.4)') niter(1) + iter, err

      Y = Y + dY * ramp2(iter, niter(2), 0.5d0)

      y_rho = merge(y_rho, 0.0_dp, ieee_is_normal(y_rho) .and. (y_rho > 0))
      y_temp = merge(y_temp, teff / 2, y_temp > teff / 2)

      if (cfg_write_all_iters) call saveiter(niter(1) + iter)

    end do relx_opacity_full

    nitert = nitert + niter(2)
  end if

  !----------------------------------------------------------------------------!

  write (ulog,*) '--- DONE'

  open(33, file = trim(outfn) // '.dat', action = 'write', status = 'new')
  call saveresult(33)
  close(33)

  !----------------------------------------------------------------------------!


  open(newunit = upar, file = trim(outfn) // '.txt', &
      action = 'write', status = 'new')
  write (upar, fmpari) "model", model
  write (upar, fmpari) "niter", nitert
  write (upar, fmpare) "zscale", fzscale(mbh,mdot,radius)
  write (upar, fmparec) "rho_0", rhoc, "Central density"
  write (upar, fmparec) "temp_0", Tc, "Central gas temperature"
  write (upar, fmparec) "Hdisk", Hdisk, "Disk height [cm]"
  write (upar, fmparfc) "alpha", alpha, "Alpha parameter"
  write (upar, fmparl) "has_corona", .FALSE.
  write (upar, fmparl) "has_magnetic", .FALSE.
  write (upar, fmparl) "has_conduction", .FALSE.
  call wpar_gl(upar)
  close(upar)

  open(34, file = trim(outfn) // ".col", action = 'write', status = 'new')
  write(34, fmcol) 'i', 'i4'
  write(34, fmcol) 'z', 'f4'
  write(34, fmcol) 'h', 'f4'
  write(34, fmcol) 'tau', 'f4'
  write(34, fmcol) 'rho', 'f4'
  write(34, fmcol) 'temp', 'f4'
  write(34, fmcol) 'frad', 'f4'
  write(34, fmcol) 'ksct', 'f4'
  write(34, fmcol) 'kabs', 'f4'
  write(34, fmcol) 'kcnd', 'f4'
  write(34, fmcol) 'pgas', 'f4'
  write(34, fmcol) 'prad', 'f4'
  close(34)

  deallocate(x,x0,Y,M,dY,tau)

contains

  !----------------------------------------------------------------------------!
  subroutine saveiter(iter)
    use relaxutils
    integer, intent(in) :: iter
    write (fn,'(A,".",I0.3,".dat")') trim(outfn),iter
    open(33, file = trim(fn), action = 'write', status = 'new')
    call saveresult(33)
    close(33)
  end subroutine

  subroutine saveresult(u)
    integer, intent(in) :: u
    integer :: i
    real(dp) :: pgas,prad
    call mktaues(x,y_rho,y_temp,tau)
    do i = 1,ngrid
      pgas =  cgs_k_over_mh / miu * y_rho(i) * y_temp(i)
      prad = cgs_a / 3 * y_trad(i)**4
      write (u,'(I6,11Es12.4)') i, x(i), x(i) / zscale, tau(i), &
      y_rho(i), y_temp(i), y_frad(i),   &
      fksct(y_rho(i), y_temp(i)), fkabs(y_rho(i), y_temp(i)), &
      fkcnd(y_rho(i), y_temp(i)), pgas, prad
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

  end subroutine

end program dv_alpha_relax
