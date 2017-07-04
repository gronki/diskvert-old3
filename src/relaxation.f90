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

    elemental subroutine mrx_sel_dims (nr, n, nbl, nbr)
        integer, intent(in) :: nr
        integer, intent(out) :: n, nbl, nbr
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

    subroutine mrx_sel_ptrs (nr, f, fbl, fbr)
        integer, intent(in) :: nr
        procedure(fcoeff_t), pointer, intent(out) :: f
        procedure(fbound_t), pointer, intent(out) :: fbl, fbr
        include 'mrxptrs.fi'
    end subroutine

    pure subroutine mrx_sel_hash (nr, ihash)
        integer, intent(in) :: nr
        integer, dimension(6), intent(out) :: ihash
        include 'mrxhash.fi'
    end subroutine

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
        integer :: i, nx,nbl,nbr,n
        procedure(fcoeff_t), pointer :: ff
        procedure(fbound_t), pointer :: fbl, fbr
        ! nieprzezroczystosci: (pochodna,rodzaj)
        ! kolejnosc: abs, sct, cond
        real(r64), dimension(3,3) :: FV
        integer, dimension(6) :: CC

        nx = size(x)
        call mrx_sel_dims(nr,n,nbl,nbr)
        allocate( ym(n), dy(n), MY(n,n), MD(n,n) )

        call mrx_sel_ptrs(nr,ff,fbl,fbr)
        call mrx_sel_hash(nr,CC)

        M = 0

        associate (xbl => x(1), YBL => Y(1:n),      &
                        & BL  => A(1:nbl),          &
                        & MBL => M(1:nbl,1:n))

            call kappabs(YBL(CC(1)), YBL(CC(2)), FV(1,1), FV(2,1), FV(3,1))
            call kappsct(YBL(CC(1)), YBL(CC(2)), FV(1,2), FV(2,2), FV(3,2))
            call kappcnd(YBL(CC(1)), YBL(CC(2)), FV(1,3), FV(2,3), FV(3,3))
            call fbl(xbl, YBL, FV, BL, MBL)

        end associate

        associate (xbr => x(nx), YBR => Y((nx-1)*n+1:), &
                        & BR => A(nbl+(nx-1)*n+1:),     &
                        & MBR => M(nbl+(nx-1)*n+1:, (nx-1)*n+1:))

            call kappabs(YBR(CC(1)), YBR(CC(2)), FV(1,1), FV(2,1), FV(3,1))
            call kappsct(YBR(CC(1)), YBR(CC(2)), FV(1,2), FV(2,2), FV(3,2))
            call kappcnd(YBR(CC(1)), YBR(CC(2)), FV(1,3), FV(2,3), FV(3,3))
            call fbr(xbr, YBR, FV, BR, MBR)

        end associate

        through_space: do i = 1, nx-1

            associate ( dx => x(i+1) - x(i), xm => (x(i+1) + x(i)) / 2, &
                    &   Ai => A(nbl+(i-1)*n+1:nbl+i*n),                 &
                    &   M1 => M(nbl+(i-1)*n+1:nbl+i*n, (i-1)*n+1:i*n),  &
                    &   M2 => M(nbl+(i-1)*n+1:nbl+i*n, i*n+1:(i+1)*n),  &
                    &   Y1 => Y((i-1)*n+1:i*n), Y2 => Y(i*n+1:(i+1)*n))

                YM = (Y2 + Y1) / 2
                DY = (Y2 - Y1) / dx

                call kappabs(YM(CC(1)), YM(CC(2)), FV(1,1), FV(2,1), FV(3,1))
                call kappsct(YM(CC(1)), YM(CC(2)), FV(1,2), FV(2,2), FV(3,2))
                call kappcnd(YM(CC(1)), YM(CC(2)), FV(1,3), FV(2,3), FV(3,3))

                call ff(xm, YM, DY, FV, Ai, MY, MD)

                M1(:,:) = MY / 2 - MD / dx
                M2(:,:) = MY / 2 + MD / dx
            end associate

        end do through_space

        deallocate(MY, MD, YM, DY)

    end subroutine

    subroutine mrx_advance(nr, x, Y, M, dY)

        integer, intent(in) :: nr
        real(r64), dimension(:), intent(in) :: x
        real(r64), dimension(:), intent(in) :: Y
        real(r64), dimension(:), intent(out), contiguous :: dY
        real(r64), dimension(:,:), intent(out), contiguous  :: M
        integer, dimension(size(dY)) :: ipiv
        integer :: errno

        call mrx_matrix(nr, x, Y, M, dY)
        call dgesv(size(M,2), 1, M, size(M,1), ipiv, dY, size(dY), errno)

        if ( errno > 0 ) then
            write (uerr, '("Parameter ", I0, " had illegal value")') abs(errno)
            error stop "zero on matrix diagonal: singular matrix"
        end if
        if ( errno < 0 ) then
            write (uerr, '("Parameter ", I0, " had illegal value")') abs(errno)
            error stop "illegal parameter value"
        end if

    end subroutine

  !----------------------------------------------------------------------------!
  ! this is a very smooth, s-shaped ramp
  ! slow convergence but good if we are away from the solution
  elemental function ramp1(i,n,y0) result(y)
    integer, intent(in) :: i,n
    real(r64), intent(in), optional :: y0
    real(r64) :: t,y
    t = merge(real(i) / n, 1.0, i .le. n)
    y = (3 - 2*t) * t**2
    if (present(y0)) y = y0 + (1 - y0) * y
  end function

  !----------------------------------------------------------------------------!
  ! this is a very steep ramp, for later iterations
  elemental function ramp2(i,n,y0) result(y)
    integer, intent(in) :: i,n
    real(r64), intent(in), optional :: y0
    real(r64) :: t,y
    t = merge(real(i) / n, 1.0, i .le. n)
    y = t * (2 - t)
    if (present(y0)) y = y0 + (1 - y0) * y
  end function

  include 'coefficients.fi'

end module relaxation
