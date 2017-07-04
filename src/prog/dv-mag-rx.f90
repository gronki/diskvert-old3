program dv_mag_relax

  use confort
  use globals
  use settings
  use iso_fortran_env, only: sp => real32, dp => real64
  use fileunits
  use relaxation
  use ss73solution, only: apxdisk, apx_estim, apx_refine
  use grid, only: space_linear, space_linlog

  !----------------------------------------------------------------------------!
  implicit none

  type(config) :: cfg
  integer :: model
  integer :: ny = 3, i, iter, nitert = 0
  integer, dimension(3) :: niter = [ 36, 12, 52 ]
  character(2**8) :: fn
  real(dp), allocatable, target :: x(:), x0(:), Y(:), dY(:), M(:,:)
  real(dp), pointer, dimension(:) :: y_rho, y_temp, y_frad, y_pmag, y_Trad
  real(dp) :: rhoc, Tc, Hdisk, err, t
  integer, dimension(6) :: c_

  !----------------------------------------------------------------------------!

  ngrid = 2**7
  htop = 18

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
  model = mrx_number( .FALSE., .TRUE., .FALSE. )
  ny = mrx_ny(model)
  call mrx_sel_hash(model, C_)

  allocate( x(ngrid), x0(ngrid), Y(ny*ngrid), dY(ny*ngrid), M(ny*ngrid,ny*ngrid) )

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
  call saveiter(0)

  kramers_opacity_ff = .false.
  kramers_opacity_bf = .false.
  relx_opacity_es : do iter = 1,niter(1)

    call mrx_advance(model, x, Y, M, dY)

    err = sum(merge((dY/Y)**2, 0d0, Y.ne.0)) / sum(merge(1, 0, Y.ne.0))
    write(ulog,'(I5,Es12.4)') nitert+iter, err

    Y = Y + dY * ramp1(iter,niter(1))
    y_rho = merge(y_rho, epsilon(1d0)*rhoc, y_rho > epsilon(1d0)*rhoc)
    y_temp = merge(y_temp, teff / 2, y_temp > teff / 2)

    call saveiter(nitert + iter)

  end do relx_opacity_es

  nitert = nitert + niter(1)

  !----------------------------------------------------------------------------!
  write (ulog,*) '--- ff+bf opacity is on'

  kramers_opacity_ff = .true.
  kramers_opacity_bf = .false.
  relx_opacity_full : do iter = 1,niter(2)

    call mrx_advance(model, x, Y, M, dY)

    err = sum(merge((dY/Y)**2, 0d0, Y.ne.0)) / sum(merge(1, 0, Y.ne.0))
    write(ulog,'(I5,Es12.4)') nitert+iter, err

    Y = Y + dY * ramp2(iter,niter(2))
    y_rho = merge(y_rho, epsilon(1d0)*rhoc, y_rho > epsilon(1d0)*rhoc)
    y_temp = merge(y_temp, teff / 2, y_temp > teff / 2)

    call saveiter(nitert + iter)

  end do relx_opacity_full

  nitert = nitert + niter(2)

  !----------------------------------------------------------------------------!
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

  relx_corona : do iter = 1,niter(3)

    call mrx_advance(model, x, Y, M, dY)

    err = sum(merge((dY/Y)**2, 0d0, Y.ne.0)) / sum(merge(1, 0, Y.ne.0))
    write(ulog,'(I5,Es12.4)') nitert+iter, err

    Y = Y + dY * ramp1(iter,niter(3)) * 0.01
    y_rho = merge(y_rho, epsilon(1d0)*rhoc, y_rho > epsilon(1d0)*rhoc)
    y_trad = merge(y_trad, teff / 2, y_trad > teff / 2)
    y_temp = merge(y_temp, teff / 2, y_temp > teff / 2)

    call saveiter(nitert + iter)

  end do relx_corona

  write (ulog,*) '--- DONE'

  open(newunit = upar, file = trim(outfn) // '.txt', action = 'write')
  write (upar, fmparec) "rho_0", rhoc, "Central density"
  write (upar, fmparec) "temp_0", Tc, "Central gas temperature"
  write (upar, fmparec) "Hdisk", Hdisk, "Disk height [cm]"
  write (upar, fmparfc) "alpha", alpha, "Alpha parameter"
  write (upar, fmparfc) "zeta", zeta, "Zeta parameter"
  call wpar_gl(upar)
  close(upar)

  deallocate(x,x0,Y,M,dY)

contains

  !----------------------------------------------------------------------------!
  subroutine mrx_transfer(nr, newnr, z, Y)
    use slf_cgs, only: cgs_k_over_mh
    integer, intent(in) :: newnr
    integer, intent(inout) :: nr
    real(dp), intent(in), dimension(:) :: z
    real(dp), intent(inout), dimension(:), allocatable :: Y
    real(dp), dimension(:), allocatable :: y_old
    ! column order: rho, Tgas, Trad, Frad, Pmag, Fcond
    integer, dimension(6) :: c_, c_old_
    integer :: ny, ny_old

    call mrx_sel_hash(nr,c_old_)
    ny_old = mrx_ny(nr)
    call mrx_sel_hash(newnr,c_)
    ny = mrx_ny(newnr)

    call move_alloc(y, y_old)
    allocate(y(size(z)*ny))

    ! density, temperature and radiative flux are present in all models
    Y(c_(1)::ny) = Y_old(c_old_(1)::ny_old)
    Y(c_(3)::ny) = Y_old(c_old_(3)::ny_old)
    Y(c_(4)::ny) = Y_old(c_old_(4)::ny_old)

    ! if radiative and gas temperature are different in the new model,
    ! then set it equal to "diffusive" temperature
    if ( c_(2) .ne. c_(3) ) then
      Y(c_(2)::ny) = Y_old(c_old_(2)::ny_old)
    end if

    ! transfer magnetic pressure if present in new model
    if ( c_(5) .ne. 0 ) then
      ! if the old model also had it, just copy
      if ( c_old_(5) .ne. 0 ) then
        Y(c_(5)::ny) = Y_old(c_old_(5)::ny_old)
      else
        ! if not, set magnetic beta to constant value of 100
        Y(c_(5)::ny) = 0.01 * 2 * cgs_k_over_mh * Y(c_(1)::ny) * Y(c_(2)::ny)
      end if
    end if

    ! transfer thermal conduction flux if present in new model
    if ( c_(6) .ne. 0 ) then
      ! if the old model also had it, just copy
      if ( c_old_(6) .ne. 0 ) then
        Y(c_(6)::ny) = Y_old(c_old_(6)::ny_old)
      else
        ! if not, just set to zero
        Y(c_(6)::ny) = 0
      end if
    end if

    deallocate(y_old)
    nr = newnr

  end subroutine

  !----------------------------------------------------------------------------!
  subroutine saveiter(iter)
    integer, intent(in) :: iter
    integer :: i
    write (fn,'(A,".",I0.3,".dat")') trim(outfn),iter
    open(33, file = trim(fn), action = 'write', status = 'new')
    do i = 1,ngrid
      write (33,'(I6,7Es12.4)') i, x(i), &
            & y_rho(i), y_temp(i), y_trad(i), y_frad(i), y_pmag(i), 0d0
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

    call mincf_get(cfg, "zeta", buf)
    if ( mincf_failed(cfg) )  then
      error stop "Magnetic field parameter (key: zeta) is REQUIRED!"
    end if
    read (buf,*) zeta

  end subroutine

end program
