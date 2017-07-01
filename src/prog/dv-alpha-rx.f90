program dv_alpha_relax

  use confort
  use globals
  use settings
  use iso_fortran_env, only: sp => real32, dp => real64
  use fileunits
  use relaxation
  use ss73solution, only: apxdisk, dvapx1, diskapx2
  use grid, only: space_linear, space_linlog

  !----------------------------------------------------------------------------!
  implicit none

  type(config) :: cfg
  integer :: model
  integer :: ny = 3, i, iter
  integer, dimension(2) :: niter = [ 18, 8 ]
  character(2**8) :: fn
  real(dp), allocatable, target :: x(:), Y(:), dY(:), M(:,:)
  real(dp), pointer, dimension(:) :: y_rho, y_temp, y_frad
  real(dp) :: rhoc, Tc, Hdisk, err

  ngrid = 2**8

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

  !----------------------------------------------------------------------------!
  allocate( x(ngrid), Y(ny*ngrid), dY(ny*ngrid), M(ny*ngrid,ny*ngrid) )

  forall (i = 1:ngrid)
    x(i) = space_linear(i,ngrid,12d0) * zscale
  end forall

  y_rho   => Y(1::ny)
  y_temp  => Y(2::ny)
  y_frad  => Y(3::ny)

  !----------------------------------------------------------------------------!
  call dvapx1(mbh, mdot, radius, alpha, rhoc, Tc, Hdisk)
  call diskapx2(mbh, mdot, radius, alpha, rhoc, Tc, Hdisk)
  call apxdisk(rhoc, Tc, Teff, Hdisk, x, y_rho, y_temp, y_frad)

  !----------------------------------------------------------------------------!
  call saveiter(0)

  kramers_opacity_ff = .false.
  kramers_opacity_bf = .false.
  relx_opacity_es : do iter = 1,niter(1)

    call mrx_advance(model, x, Y, M, dY)

    err = sum(merge((dY/Y)**2, 0d0, Y.ne.0)) / sum(merge(1, 0, Y.ne.0))
    write(ulog,'(I5,Es12.4)') iter, err

    Y = Y + dY * ramp1(iter,niter(1))
    y_rho = merge(y_rho, 0.0_dp, y_rho > 0)
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
    y_rho = merge(y_rho, 0.0_dp, y_rho > 0)
    y_temp = merge(y_temp, teff / 2, y_temp > teff / 2)

    call saveiter(niter(1) + iter)

  end do relx_opacity_full
  !----------------------------------------------------------------------------!

  write (ulog,*) '--- DONE'

  deallocate(x,Y,M,dY)

contains

  !----------------------------------------------------------------------------!
  subroutine saveiter(iter)
    integer, intent(in) :: iter
    integer :: i
    write (fn,'(A,".",I0.3,".dat")') trim(outfn),iter
    open(33, file = trim(fn), action = 'write', status = 'new')
    do i = 1,ngrid
      write (33,'(I6,4Es12.4)') i, x(i) / zscale, y_rho(i) / rhoc, y_temp(i) / teff, y_frad(i) / facc
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
