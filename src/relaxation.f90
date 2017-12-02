module relaxation

  use iso_fortran_env, only: r64 => real64
  use ieee_arithmetic
  use fileunits

  use slf_cgs
  use globals

  implicit none

  interface
    pure subroutine funout_t(z, Y, YY)
      import r64
      real(r64), intent(in) :: z
      real(r64), intent(in), dimension(:) :: Y
      real(r64), intent(out), dimension(:) :: YY
    end subroutine
    pure subroutine fun0_t(z, Y, F, A, M)
      import r64
      real(r64), intent(in) :: z
      real(r64), intent(in), dimension(:) :: Y
      real(r64), intent(in), dimension(:,:) :: F
      real(r64), intent(out), dimension(:) :: A
      real(r64), intent(out), dimension(:,:) :: M
    end subroutine
    pure subroutine fun1_t(z, Y, D, F, A, MY, MD)
      import r64
      real(r64), intent(in) :: z
      real(r64), intent(in), dimension(:) :: Y,D
      real(r64), intent(in), dimension(:,:) :: F
      real(r64), intent(out), dimension(:) :: A
      real(r64), intent(out), dimension(:,:) :: MY,MD
    end subroutine
    pure subroutine fun2_t(z, Y, D, D2, F, A, MY, MD, MD2)
      import r64
      real(r64), intent(in) :: z
      real(r64), intent(in), dimension(:) :: Y,D,D2
      real(r64), intent(in), dimension(:,:) :: F
      real(r64), intent(out), dimension(:) :: A
      real(r64), intent(out), dimension(:,:) :: MY,MD,MD2
    end subroutine
  end interface

  real(r64) :: alpha = 0, zeta = 0, nu = 0
  real(r64) :: omega, radius, facc, teff, zscale

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

  pure function mrx_number(bil, magnetic, conduction) result(nr)
    logical, intent(in) :: magnetic,conduction
    character, intent(in) :: bil
    integer :: nr

    nr = 1

    select case (bil)
    case (EQUATION_DIFFUSION)
    case (EQUATION_COMPTON)
      nr = nr + 1
    case (EQUATION_BALANCE)
      nr = nr + 3
    case default
      error stop "wrong character for bil, allowed: D,W,C"
    end select

    if (magnetic)   nr = nr + 4
    if (conduction) nr = nr + 8

  end function

  elemental subroutine mrx_sel_dims (nr, ny, neq0, neq1, nbl, nbr)
    integer, intent(in) :: nr
    integer, intent(out) :: ny, neq1, nbl, nbr, neq0
    include 'mrxdims.fi'
  end subroutine

  elemental subroutine mrx_sel_nvar (nr, ny)
    integer, intent(in) :: nr
    integer, intent(out) :: ny
    integer :: neq1, nbl, nbr, neq0
    call mrx_sel_dims (nr, ny, neq0, neq1, nbl, nbr)
  end subroutine

  subroutine mrx_sel_funcs (nr, feq0, feq1, fbl, fbr, fout)
    integer, intent(in) :: nr
    procedure(fun1_t), pointer, intent(out) :: feq1
    procedure(funout_t), pointer, intent(out) :: fout
    procedure(fun0_t), pointer, intent(out) :: fbl, fbr, feq0
    include 'mrxptrs.fi'
  end subroutine

  pure subroutine mrx_sel_hash (nr, ihash)
    integer, intent(in) :: nr
    integer, dimension(:), intent(out) :: ihash
    include 'mrxhash.fi'
  end subroutine

  pure subroutine mrx_sel_fout (nr, fout)
    integer, intent(in) :: nr
    procedure(funout_t), pointer, intent(out) :: fout
    include 'mrxfout.fi'
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
    integer :: i, nx,nbl,nbr,ny,neq1,neq0
    procedure(fun1_t), pointer :: feq1
    procedure(fun0_t), pointer :: fbl, fbr, feq0
    procedure(funout_t), pointer :: fout
    ! nieprzezroczystosci: (pochodna,rodzaj)
    ! kolejnosc: abs, sct, cond
    real(r64), dimension(3,3) :: FV
    integer, dimension(6) :: c_

    nx = size(x)
    call mrx_sel_dims(nr,ny,neq0,neq1,nbl,nbr)
    allocate( ym(ny), dy(ny), MY(neq1,ny), MD(neq1,ny) )

    call mrx_sel_funcs(nr,feq0,feq1,fbl,fbr,fout)
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
                & BR  => A(nbl+(nx-1)*ny+1+neq0:nx*ny),           &
                & MBR => M(nbl+(nx-1)*ny+1+neq0:nx*ny, (nx-1)*ny+1:nx*ny))

        call kappabs(YBR(c_(1)), YBR(c_(2)), FV(1,1), FV(2,1), FV(3,1))
        call kappsct(YBR(c_(1)), YBR(c_(2)), FV(1,2), FV(2,2), FV(3,2))
        call kappcnd(YBR(c_(1)), YBR(c_(2)), FV(1,3), FV(2,3), FV(3,3))

        call FBR(xbr, YBR, FV, BR, MBR)

      end associate
    end if

    if ( neq0 > 0 ) then
      fillc: do i = 1,nx

        associate ( xc => x(i), YC => Y((i-1)*ny+1:i*ny),    &
                  &  C => A(nbl+(i-1)*ny+1:nbl+(i-1)*ny+neq0), &
                  & MC => M(nbl+(i-1)*ny+1:nbl+(i-1)*ny+neq0, (i-1)*ny+1:i*ny))

          call kappabs(YC(c_(1)), YC(c_(2)), FV(1,1), FV(2,1), FV(3,1))
          call kappsct(YC(c_(1)), YC(c_(2)), FV(1,2), FV(2,2), FV(3,2))
          call kappcnd(YC(c_(1)), YC(c_(2)), FV(1,3), FV(2,3), FV(3,3))

          call feq0(xc, YC, FV, C, MC)

        end associate
      end do fillc
    end if

    filla: do i = 1, nx-1

      associate ( dx => x(i+1) - x(i), xm => (x(i+1) + x(i)) / 2, &
              &   Ai => A(nbl+(i-1)*ny+1+neq0:nbl+i*ny),                 &
              &   M1 => M(nbl+(i-1)*ny+1+neq0:nbl+i*ny, (i-1)*ny+1:i*ny),  &
              &   M2 => M(nbl+(i-1)*ny+1+neq0:nbl+i*ny, i*ny+1:(i+1)*ny),  &
              &   Y1 => Y((i-1)*ny+1:i*ny), Y2 => Y(i*ny+1:(i+1)*ny))

        YM(:) = (Y2 + Y1) / 2
        DY(:) = (Y2 - Y1) / dx

        call kappabs(YM(c_(1)), YM(c_(2)), FV(1,1), FV(2,1), FV(3,1))
        call kappsct(YM(c_(1)), YM(c_(2)), FV(1,2), FV(2,2), FV(3,2))
        call kappcnd(YM(c_(1)), YM(c_(2)), FV(1,3), FV(2,3), FV(3,3))

        call feq1(xm, YM, DY, FV, Ai, MY, MD)

        M1(:,:) = MY / 2 - MD / dx
        M2(:,:) = MY / 2 + MD / dx

      end associate

    end do filla

    deallocate(MY, MD, YM, DY)

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

  ! this is a very smooth ramp where you can control the starting
  ! level as well as slope of the starting point. It may be applied
  ! for very unstable equations (A = 0, B = 0) or for a very quick
  ! convergence (A = 0.25, B = 1)
  elemental function ramp4(i,n,A,B) result(y)
    integer, intent(in) :: i,n
    real(r64), intent(in) :: A,B
    real(r64) :: x,y
    x = merge(real(i) / n, 1.0, i .le. n)
    y = A + (1 - A) * B * x           &
      + 3 * (A - 1) * (B - 2) * x**2  &
      - (A - 1) * (3*B - 8) * x**3    &
      + (A - 1) * (B - 3) * x**4
  end function

  !----------------------------------------------------------------------------!

  subroutine mrx_relax(nr, x, Y, niter, r, errno)

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

    errno = 0
    call mrx_sel_hash(nr,ch)
    call mrx_sel_nvar(nr,ny)
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
  subroutine mrx_transfer(nr_old, nr, z, Y)
    use slf_cgs, only: cgs_k_over_mh
    integer, intent(in) :: nr
    integer, intent(inout) :: nr_old
    real(r64), intent(in), dimension(:) :: z
    real(r64), intent(inout), dimension(:), allocatable, target :: Y
    real(r64), dimension(:), allocatable, target :: y_old
    real(r64), dimension(:,:), pointer :: yv, yv_old
    ! column order: rho, Tgas, Trad, Frad, Pmag, Fcond
    integer, dimension(6) :: c_, c_old_
    integer :: ny, ny_old

    call mrx_sel_hash(nr_old, c_old_)
    call mrx_sel_nvar(nr_old, ny_old)

    call mrx_sel_hash(nr, c_)
    call mrx_sel_nvar(nr, ny)

    call move_alloc(y, y_old)
    allocate(y(size(z)*ny))

    yv(1:ny,1:size(z)) => y
    yv_old(1:ny_old,1:size(z)) => y_old

    ! density, temperature and radiative flux are present in all models
    yv(c_(1),:) = yv_old(c_old_(1),:)
    yv(c_(2),:) = yv_old(c_old_(2),:)
    yv(c_(4),:) = yv_old(c_old_(4),:)

    ! if radiative and gas temperature are different in the new model,
    ! then set it equal to "diffusive" temperature
    if ( c_(2) .ne. c_(3) ) then
      yv(c_(3),:) = yv_old(c_old_(3),:)
    end if

    ! transfer magnetic pressure if present in new model
    if ( c_(5) .ne. 0 ) then
      ! if the old model also had it, just copy
      if ( c_old_(5) .ne. 0 ) then
        yv(c_(5),:) = yv_old(c_old_(5),:)
      else
        ! if not, set magnetic beta to constant value of 100
        yv(c_(5),:) = 0.01 * 2 * cgs_k_over_mh * yv(c_(1),:) * yv(c_(2),:)
      end if
    end if

    ! transfer thermal conduction flux if present in new model
    if ( c_(6) .ne. 0 ) then
      ! if the old model also had it, just copy
      if ( c_old_(6) .ne. 0 ) then
        yv(c_(6),:) = yv_old(c_old_(6),:)
      else
        ! if not, just set to zero
        yv(c_(6),:) = 0
      end if
    end if

    deallocate(y_old)
    nr_old = nr

  end subroutine

  !----------------------------------------------------------------------------!

  include 'mrxcoeff.fi'

end module relaxation
