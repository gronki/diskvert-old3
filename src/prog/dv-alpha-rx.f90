program dv_alpha_relax

  use confort
  use globals
  use settings
  use iso_fortran_env, only: sp => real32, dp => real64
  use ieee_arithmetic, only: ieee_is_normal, ieee_is_nan
  use fileunits
  use relaxation
  use ss73solution, only: apxdisk, apx_estim, apx_refine
  use grid, only: space_linear, space_linlog

  !----------------------------------------------------------------------------!
  implicit none

  type(config) :: cfg
  integer :: model
  integer :: ny = 3, i, iter
  integer, dimension(2) :: niter = [ 36, 8 ]
  character(2**8) :: fn
  real(dp), allocatable, target :: x(:), x0(:), Y(:), dY(:), M(:,:)
  real(dp), pointer, dimension(:) :: y_rho, y_temp, y_frad
  real(dp) :: rhoc, Tc, Hdisk, err, t
  integer, dimension(6) :: C_


  !----------------------------------------------------------------------------!

  ngrid = 2**7
  htop = 4

  !----------------------------------------------------------------------------!

  call rdargvgl
  call mincf_read(cfg)
  call rdconfgl(cfg)
  call rdconf(cfg)
  call mincf_free(cfg)

  !----------------------------------------------------------------------------!
  ! calculate the global parameters
  call mrx_init(mbh, mdot, radius, alpha)

  ! get the model number
  model = mrx_number( .FALSE., .FALSE., .FALSE. )
  ny = mrx_ny(model)
  call mrx_sel_hash(model, C_)

  allocate( x(ngrid), x0(ngrid), Y(ny*ngrid), dY(ny*ngrid), M(ny*ngrid,ny*ngrid) )

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
  y_frad  => Y(C_(4)::ny)

  !----------------------------------------------------------------------------!

  y_frad = x0 * facc
  y_temp = exp(-0.5*(x/Hdisk)**2) * (Tc - 0.841 * Teff) + 0.841 * Teff
  y_rho =  rhoc * (exp(-0.5*(x/Hdisk)**2) + 1e-5)

  !----------------------------------------------------------------------------!
  call saveiter(0)

  kramers_opacity_ff = .false.
  kramers_opacity_bf = .false.
  relx_opacity_es : do iter = 1,niter(1)

    call mrx_advance(model, x, Y, M, dY)

    err = sum(merge((dY/Y)**2, 0d0, Y.ne.0)) / sum(merge(1, 0, Y.ne.0))
    write(ulog,'(I5,Es12.4)') iter, err

    Y = Y + dY * ramp1(iter,niter(1))

    y_rho = merge(y_rho, 0.0_dp, ieee_is_normal(y_rho) .and. (y_rho > 0))
    y_temp = merge(y_temp, teff / 2, y_temp > teff / 2)

    call saveiter(iter)

  end do relx_opacity_es

  !----------------------------------------------------------------------------!
  write (ulog,*) '--- bf+ff opacity is ON'

  kramers_opacity_ff = .true.
  kramers_opacity_bf = .true.
  relx_opacity_full : do iter = 1,niter(2)

    call mrx_advance(model, x, Y, M, dY)

    err = sum(merge((dY/Y)**2, 0d0, Y.ne.0)) / sum(merge(1, 0, Y.ne.0))
    write(ulog,'(I5,Es12.4)') niter(1) + iter, err

    Y = Y + dY * ramp2(iter, niter(2), 0.5d0)

    y_rho = merge(y_rho, 0.0_dp, ieee_is_normal(y_rho) .and. (y_rho > 0))
    y_temp = merge(y_temp, teff / 2, y_temp > teff / 2)

    call saveiter(niter(1) + iter)

  end do relx_opacity_full
  !----------------------------------------------------------------------------!

  write (ulog,*) '--- DONE'

  open(newunit = upar, file = trim(outfn) // '.txt', action = 'write')
  write (upar, fmparec) "rho_0", rhoc, "Central density"
  write (upar, fmparec) "temp_0", Tc, "Central gas temperature"
  write (upar, fmparec) "Hdisk", Hdisk, "Disk height [cm]"
  write (upar, fmparfc) "alpha", alpha, "Alpha parameter"
  call wpar_gl(upar)
  close(upar)

  deallocate(x,x0,Y,M,dY)

contains

  !----------------------------------------------------------------------------!
  subroutine saveiter(iter)
    integer, intent(in) :: iter
    integer :: i
    write (fn,'(A,".",I0.3,".dat")') trim(outfn),iter
    open(33, file = trim(fn), action = 'write', status = 'new')
    do i = 1,ngrid
      write (33,'(I6,7Es12.4)') i, x(i), y_rho(i), y_temp(i), y_temp(i),  &
            y_frad(i), 0d0, 0d0, fksct(y_rho(i), y_temp(i)),  &
            fkabs(y_rho(i), y_temp(i))
    end do
    close(33)
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
