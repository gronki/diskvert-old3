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

    private :: ramp

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

    pure subroutine mrx_sel_dims (nr, n, nbl, nbr)
        integer, intent(in) :: nr
        integer, intent(out) :: n, nbl, nbr
        include 'mrxdims.fi'
    end subroutine

    elemental function mrx_ny(nr) result(ny)
        integer, intent(in) :: nr
        integer :: ny
        include 'mrxynum.fi'
    end function

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
        real(r64), dimension(3,3) :: fval

        nx = size(x)
        call mrx_sel_dims(nr,n,nbl,nbr)
        allocate( ym(n), dy(n), MY(n,n), MD(n,n) )

        call mrx_sel_ptrs(nr,ff,fbl,fbr)

        M = 0

        associate (xbl => x(1), YBL => Y(1:n),      &
                        & BL  => A(1:nbl),          &
                        & MBL => M(1:nbl,1:n))

            call fbl(xbl, YBL, fval, BL, MBL)

        end associate

        associate (xbr => x(nx), YBR => Y((nx-1)*n+1:), &
                        & BR => A(nbl+(nx-1)*n+1:),     &
                        & MBR => M(nbl+(nx-1)*n+1:, (nx-1)*n+1:))

            call fbr(xbr, YBR, fval, BR, MBR)

        end associate

        through_space: do i = 1, nx-1

            associate ( dx => x(i+1) - x(i), xm => (x(i+1) + x(i)) / 2, &
                    &   Ai => A(nbl+(i-1)*n+1:nbl+i*n),                 &
                    &   M1 => M(nbl+(i-1)*n+1:nbl+i*n, (i-1)*n+1:i*n),  &
                    &   M2 => M(nbl+(i-1)*n+1:nbl+i*n, i*n+1:(i+1)*n),  &
                    &   Y1 => Y((i-1)*n+1:i*n), Y2 => Y(i*n+1:(i+1)*n))

                YM = (Y2 + Y1) / 2
                DY = (Y2 - Y1) / dx
                call ff(xm, YM, DY, fval, Ai, MY, MD)

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

    elemental function ramp(i,n) result(y)
        integer, intent(in) :: i,n
        real(r64) :: t,y
        t = merge(real(i) / n, 1.0, i .le. n)
        y = 3 * t**2 - 2 * t**3
    end function

    include 'coefficients.fi'

end module relaxation
