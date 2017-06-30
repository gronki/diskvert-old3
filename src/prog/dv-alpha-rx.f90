program dv_alpha_relax

  use confort
  use globals
  use settings
  use iso_fortran_env, only: sp => real32, dp => real64
  use fileunits
  use relaxation
  use ss73solution, only: apxdisk, dvapx1, diskapx2
  use grid, only: space_linear, space_linlog

  implicit none

  type(config) :: cfg
  integer :: model
  integer :: ny = 3, i, iter
  character(2**8) :: fn
  real(dp), allocatable, target :: x(:), Y(:), dY(:), M(:,:)
  real(dp), pointer, dimension(:) :: y_rho, y_temp, y_frad
  real(dp) :: rhoc, Tc, Hdisk

  ngrid = 2**8

  call rdargvgl
  call mincf_read(cfg)
  call rdconfgl(cfg)
  call rdconf(cfg)
  call mincf_free(cfg)

  ! calculate the global parameters
  call mrx_init(mbh, mdot, radius, alpha)

  ! get the model number
  model = mrx_number( .FALSE., .FALSE., .FALSE. )
  ny = mrx_ny(model)

  allocate( x(ngrid), Y(ny*ngrid), dY(ny*ngrid), M(ny*ngrid,ny*ngrid) )

  forall (i = 1:ngrid)
    x(i) = space_linlog(i,ngrid,24d0) * zscale
  end forall

  y_rho   => Y(1::ny)
  y_temp  => Y(2::ny)
  y_frad  => Y(3::ny)

  call dvapx1(mbh, mdot, radius, alpha, rhoc, Tc, Hdisk)
  call diskapx2(mbh, mdot, radius, alpha, rhoc, Tc, Hdisk)
  call apxdisk(rhoc, Tc, Teff, Hdisk, x, y_rho, y_temp, y_frad)

  call saveiter(0)
  relaxate : do iter = 1,4
    call mrx_advance(model, x, Y, M, dY)
    Y = Y + dY
    call saveiter(iter)
  end do relaxate

  deallocate(x,Y,M,dY)

contains

  subroutine saveiter(iter)
    integer, intent(in) :: iter
    integer :: i
    write (fn,'(A,".",I0.3,".dat")') trim(outfn),iter
    open(33, file = trim(fn), action = 'write')
    do i = 1,ngrid
      write (33,'(I6,4Es12.4)') i, x(i), y_rho(i) / rhoc, y_temp(i) / teff, y_frad(i) / facc
    end do
    close(33)
  end subroutine

  elemental function ramp(i,n) result(y)
    integer, intent(in) :: i,n
    real(r64) :: t,y
    t = merge(real(i) / n, 1.0, i .le. n)
    y = (3 - 2*t) * t**2
  end function

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
