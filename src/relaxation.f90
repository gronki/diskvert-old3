module relaxation

  use iso_fortran_env, only: r64 => real64
  use ieee_arithmetic
  use fileunits

  use slf_cgs
  use globals

  implicit none

  abstract interface
    pure subroutine fcoeff_t(z, Y, D, F, A, MY, MD)
      import r64
      real(r64), intent(in) :: z
      real(r64), intent(in), dimension(:) :: Y,D
      real(r64), intent(in), dimension(:,:) :: F
      real(r64), intent(out), dimension(:) :: A
      real(r64), intent(out), dimension(:,:) :: MY,MD
    end subroutine
    pure subroutine fbound_t(z, Y, F, B, M)
      import r64
      real(r64), intent(in) :: z
      real(r64), intent(in), dimension(:) :: Y
      real(r64), intent(in), dimension(:,:) :: F
      real(r64), intent(out), dimension(:) :: B
      real(r64), intent(out), dimension(:,:) :: M
    end subroutine
  end interface

  real(r64) :: alpha = 0, zeta = 0
  real(r64) :: omega, radius, facc, teff, zscale

  integer, parameter :: c_rho = 1, c_Tgas = 2, c_Trad = 3, &
  & c_Frad = 4, c_Pmag = 5, c_Fcond = 6

contains

  subroutine mrx_init(mb, md, rad, alph, zet)
    real(r64), intent(in) :: mb, md, rad, alph
    real(r64), intent(in), optional :: zet
    mbh = mb
    mdot = md
    radius = rad
    alpha = alph
    if (present(zet)) zeta = zet
    call cylinder(mbh, mdot, radius, omega, facc, teff, zscale)
  end subroutine

  !----------------------------------------------------------------------------!

  pure function mrx_number(compton, magnetic, conduction) result(nr)
    logical, intent(in) :: compton,magnetic,conduction
    integer :: nr
    nr = 1  + merge(1, 0, compton)      &
            + merge(2, 0, magnetic)     &
            + merge(4, 0, conduction)
  end function

  pure subroutine mrx_sel_name (nr, name)
    integer, intent(in) :: nr
    character(*), intent(out) :: name
    include 'mrxname.fi'
  end subroutine

  elemental subroutine mrx_sel_dims (nr, ny, na, nc, nbl, nbr)
    integer, intent(in) :: nr
    integer, intent(out) :: ny, na, nbl, nbr, nc
    include 'mrxdims.fi'
  end subroutine

  elemental function mrx_ny(nr) result(ny)
    integer, intent(in) :: nr
    integer :: ny
    include 'mrxynum.fi'
  end function

  elemental subroutine mrx_sel_ny(nr,ny)
    integer, intent(in) :: nr
    integer, intent(out) :: ny
    include 'mrxynum.fi'
  end subroutine

  subroutine mrx_sel_ptrs (nr, fa, fc, fbl, fbr)
    integer, intent(in) :: nr
    procedure(fcoeff_t), pointer, intent(out) :: fa
    procedure(fbound_t), pointer, intent(out) :: fbl, fbr, fc
    include 'mrxptrs.fi'
  end subroutine

  pure subroutine mrx_sel_hash (nr, ihash)
    integer, intent(in) :: nr
    integer, dimension(6), intent(out) :: ihash
    include 'mrxhash.fi'
  end subroutine

  !----------------------------------------------------------------------------!

  subroutine mrx_matrix(nr,x,Y,M,A)

    integer :: nr
    real(r64), dimension(:), intent(in) :: x
    ! uklad: [ Y1(1) Y2(1) Y3(1) Y1(2) Y2(2) Y3(2) ... ]
    real(r64), dimension(:), intent(in) :: Y
    ! uklad: [ A1(1) A2(1) A3(1) A1(2) A2(2) A3(2) ... ]
    real(r64), dimension(:), intent(out) :: A
    real(r64), dimension(:,:), target, intent(out) :: M

    real(r64), dimension(:,:), allocatable :: MY, MD
    real(r64), dimension(:), allocatable :: YM,DY
    integer :: i, nx,nbl,nbr,ny,na,nc
    procedure(fcoeff_t), pointer :: fa
    procedure(fbound_t), pointer :: fbl, fbr, fc
    ! nieprzezroczystosci: (pochodna,rodzaj)
    ! kolejnosc: abs, sct, cond
    real(r64), dimension(3,3) :: FV
    integer, dimension(6) :: c_

    nx = size(x)
    call mrx_sel_dims(nr,ny,na,nc,nbl,nbr)
    allocate( ym(ny), dy(ny), MY(na,ny), MD(na,ny) )

    call mrx_sel_ptrs(nr,fa,fc,fbl,fbr)
    call mrx_sel_hash(nr,c_)

    M = 0

    if (nbl > 0) then
      associate (xbl => x(1), YBL => Y(1:ny),     &
              &  BL  => A(1:nbl),                 &
              &  MBL => M(1:nbl,1:ny))

        call kappabs(YBL(c_(1)), YBL(c_(2)), FV(1,1), FV(2,1), FV(3,1))
        call kappsct(YBL(c_(1)), YBL(c_(2)), FV(1,2), FV(2,2), FV(3,2))
        call kappcnd(YBL(c_(1)), YBL(c_(2)), FV(1,3), FV(2,3), FV(3,3))

        call FBL(xbl, YBL, FV, BL, MBL)

      end associate
    end if

    if ( nbr > 0 ) then
      associate ( xbr => x(nx), YBR => Y((nx-1)*ny+1:nx*ny), &
                & BR  => A(nbl+(nx-1)*ny+1+nc:nx*ny),           &
                & MBR => M(nbl+(nx-1)*ny+1+nc:nx*ny, (nx-1)*ny+1:nx*ny))

        call kappabs(YBR(c_(1)), YBR(c_(2)), FV(1,1), FV(2,1), FV(3,1))
        call kappsct(YBR(c_(1)), YBR(c_(2)), FV(1,2), FV(2,2), FV(3,2))
        call kappcnd(YBR(c_(1)), YBR(c_(2)), FV(1,3), FV(2,3), FV(3,3))

        call FBR(xbr, YBR, FV, BR, MBR)

      end associate
    end if

    if ( nc > 0 ) then
      fillc: do i = 1,nx

        associate ( xc => x(i), YC => Y((i-1)*ny+1:i*ny),    &
                  &  C => A(nbl+(i-1)*ny+1:nbl+(i-1)*ny+nc), &
                  & MC => M(nbl+(i-1)*ny+1:nbl+(i-1)*ny+nc, (i-1)*ny+1:i*ny))

          call kappabs(YC(c_(1)), YC(c_(2)), FV(1,1), FV(2,1), FV(3,1))
          call kappsct(YC(c_(1)), YC(c_(2)), FV(1,2), FV(2,2), FV(3,2))
          call kappcnd(YC(c_(1)), YC(c_(2)), FV(1,3), FV(2,3), FV(3,3))

          call FC(xc, YC, FV, C, MC)

        end associate
      end do fillc
    end if

    filla: do i = 1, nx-1

      associate ( dx => x(i+1) - x(i), xm => (x(i+1) + x(i)) / 2, &
              &   Ai => A(nbl+(i-1)*ny+1+nc:nbl+i*ny),                 &
              &   M1 => M(nbl+(i-1)*ny+1+nc:nbl+i*ny, (i-1)*ny+1:i*ny),  &
              &   M2 => M(nbl+(i-1)*ny+1+nc:nbl+i*ny, i*ny+1:(i+1)*ny),  &
              &   Y1 => Y((i-1)*ny+1:i*ny), Y2 => Y(i*ny+1:(i+1)*ny))

        YM(:) = (Y2 + Y1) / 2
        DY(:) = (Y2 - Y1) / dx

        call kappabs(YM(c_(1)), YM(c_(2)), FV(1,1), FV(2,1), FV(3,1))
        call kappsct(YM(c_(1)), YM(c_(2)), FV(1,2), FV(2,2), FV(3,2))
        call kappcnd(YM(c_(1)), YM(c_(2)), FV(1,3), FV(2,3), FV(3,3))

        call FA(xm, YM, DY, FV, Ai, MY, MD)

        M1(:,:) = MY / 2 - MD / dx
        M2(:,:) = MY / 2 + MD / dx

      end associate

    end do filla

    deallocate(MY, MD, YM, DY)

  end subroutine

  subroutine mrx_advance(nr, x, Y, M, dY, errno)

    integer, intent(in) :: nr
    real(r64), dimension(:), intent(in) :: x
    real(r64), dimension(:), intent(in) :: Y
    real(r64), dimension(:), intent(out), contiguous :: dY
    real(r64), dimension(:,:), intent(out), contiguous  :: M
    integer, intent(out) :: errno
    integer, dimension(size(dY)) :: ipiv

    call mrx_matrix(nr, x, Y, M, dY)
    call dgesv(size(M,2), 1, M, size(M,1), ipiv, dY, size(dY), errno)

    if ( errno > 0 ) then
      write (uerr, '("Parameter ", I0, " had illegal value")') abs(errno)
      write (uerr, '("zero on matrix diagonal: singular matrix")')
    end if
    if ( errno < 0 ) then
      write (uerr, '("Parameter ", I0, " had illegal value")') abs(errno)
      write (uerr, '("illegal parameter value")')
    end if

  end subroutine

  !----------------------------------------------------------------------------!

  ! this is a linear ramp
  elemental function ramp1(i,n,y0) result(y)
    integer, intent(in) :: i,n
    real(r64), intent(in), optional :: y0
    real(r64) :: y
    y = merge(real(i) / n, 1.0, i .le. n)
    if (present(y0)) y = y0 + (1 - y0) * y
  end function

  ! this is a very steep ramp, for later iterations
  elemental function ramp2(i,n,y0) result(y)
    integer, intent(in) :: i,n
    real(r64), intent(in), optional :: y0
    real(r64) :: t,y
    t = merge(real(i) / n, 1.0, i .le. n)
    y = t * (2 - t)
    if (present(y0)) y = y0 + (1 - y0) * y
  end function

  ! this is a very smooth, s-shaped ramp
  ! slow convergence but good if we are away from the solution
  elemental function ramp3(i,n,y0) result(y)
    integer, intent(in) :: i,n
    real(r64), intent(in), optional :: y0
    real(r64) :: t,y
    t = merge(real(i) / n, 1.0, i .le. n)
    y = (3 - 2*t) * t**2
    if (present(y0)) y = y0 + (1 - y0) * y
  end function

  !----------------------------------------------------------------------------!

  subroutine mrx_relax(nr, x, Y, niter, r, errno)

    ! interface callback
    !   subroutine callback
    !     integer, intent(in)
    !   end subroutine callback
    ! end interface callback

    integer, intent(in) :: nr, niter
    real(r64), intent(in), dimension(:) :: x
    real(r64), intent(inout), dimension(:) :: Y
    real(r64), dimension(size(Y),size(Y)) :: M
    real(r64), dimension(size(Y)) :: dY, Yerr
    real(r64), intent(in) :: r
    real(r64) :: err
    integer, dimension(size(Y)) :: yerrmask
    integer, intent(out) :: errno
    integer :: i,it, ch(6), ny
    integer, dimension(size(Y)) :: ipiv
    ! procedure(callback) :: cb

    errno = 0
    call mrx_sel_hash(nr,ch)
    call mrx_sel_ny(nr,ny)
    write(uerr,'(A5,1X,A9)') 'ITER', 'ERROR'

    relax: do it = 1, niter

      call mrx_matrix(nr, x, Y, M, dY)
      call dgesv(size(M,2), 1, M, size(M,1), ipiv, dY, size(dY), errno)

      if (errno .ne. 0) exit relax

      where (Y .ne. 0 .and. ieee_is_normal(dy))
        yerr = (dy / y)**2
        yerrmask = 1
      elsewhere
        yerr = 0
        yerrmask = 0
      end where

      where (.not.ieee_is_normal(dy)) dy = 0

      err = sum(yerr) / sum(yerrmask)
      write(uerr,'(I5,1X,Es9.2)') it, err

      Y = Y + dY * ramp3(it,niter,r)

      if ( err .lt. 3e-8 ) then
        write (uerr, '("error = ",ES9.2,", exiting")') err
        exit relax
      end if

    end do relax

    write (uerr, '("errno: ",I0)') errno
  end subroutine

  !----------------------------------------------------------------------------!
  subroutine mrx_transfer(nr, newnr, z, Y)
    use slf_cgs, only: cgs_k_over_mh
    integer, intent(in) :: newnr
    integer, intent(inout) :: nr
    real(r64), intent(in), dimension(:) :: z
    real(r64), intent(inout), dimension(:), allocatable :: Y
    real(r64), dimension(:), allocatable :: y_old
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

  include 'coefficients.fi'

end module relaxation
