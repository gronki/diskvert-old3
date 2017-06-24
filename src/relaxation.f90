module relaxation

    use iso_fortran_env
    use ieee_arithmetic

    use precision
    use slf_cgs
    use globals
    use lapack

    ! wspolczynniki wygenerowane
    use relax_coefficients

    implicit none

    abstract interface
        pure subroutine coeff_driver_t(z, Y, D, A, MY, MD, ny)
            import fp
            integer, intent(in) :: ny
            real(fp), intent(in) :: z
            real(fp), intent(in), dimension(ny) :: Y,D
            real(fp), intent(out), dimension(ny) :: A
            real(fp), intent(out), dimension(ny,ny) :: MY,MD
        end subroutine
        pure subroutine bound_driver_t(z, Y, B, M, ny, nb)
            import fp
            integer, intent(in) :: ny,nb
            real(fp), intent(in) :: z
            real(fp), intent(in), dimension(ny) :: Y
            real(fp), intent(out), dimension(nb) :: B
            real(fp), intent(out), dimension(nb,ny) :: M
        end subroutine
    end interface

    type :: model_t
        integer :: nr = 0
        procedure(coeff_driver_t), pointer, nopass, private :: coeff
        procedure(bound_driver_t), pointer, nopass, private :: coeff_L, coeff_R
        integer, private :: ny, nbl, nbr
        integer, dimension(6), private :: ix
        real(fp), dimension(:), pointer :: x => NULL()
        real(fp), dimension(:), pointer, contiguous  :: Y,dY
        real(fp), dimension(:,:), pointer, contiguous :: M
    contains
        procedure :: init    => model_initialize
        procedure :: matrix  => model_matrix
        procedure :: alloc  => model_allocate
        procedure :: advance => model_advance
    end type

contains

    subroutine model_initialize (model, compton, magnetic, conduction)
        class(model_t) :: model
        logical, intent(in) :: compton,magnetic,conduction

        model % nr = 1  + merge(1, 0, compton)      &
                        + merge(2, 0, magnetic)     &
                        + merge(4, 0, conduction)

        include 'model_select.inc'

    end subroutine

    subroutine model_allocate(model,x)

        class(model_t) :: model
        real(fp), dimension(:), intent(in), target :: x

        model % x => x

        ! nl => ny + nbl - 1
        ! nu => 2 * ny - nbl - 1

        associate ( ny => model % ny,       &
                    nbl => model % nbl,     &
                    nbr => model % nbr,     &
                    nx => size(x) )

            allocate( model % Y(nx*ny), model % dY(nx*ny) )
            allocate( model % M(nx*ny,nx*ny) )

        end associate


    end subroutine


    subroutine model_matrix(model,x,Y,M,A)

        class(model_t) :: model
        real(fp), dimension(:), intent(in) :: x
        ! uklad: [ Y1(1) Y2(1) Y3(1) Y1(2) Y2(2) Y3(2) ... ]
        real(fp), dimension(size(x) * (model % ny)), intent(in) :: Y
        ! uklad: [ A1(1) A2(1) A3(1) A1(2) A2(2) A3(2) ... ]
        real(fp), dimension(size(x) * (model % ny)), intent(out) :: A
        real(fp), dimension(size(A), size(Y)), target, intent(out) :: M

        real(fp), dimension(model % ny) :: Ym, Dm
        real(fp), dimension(model % ny) :: Am
        real(fp), dimension(model % ny, model % ny) :: MY, MD
        real(fp), dimension(model % nbl) :: BL
        real(fp), dimension(model % nbl, model % ny) :: MBL
        real(fp), dimension(model % nbr) :: BR
        real(fp), dimension(model % nbr, model % ny) :: MBR

        real(fp) :: dx, xm
        integer :: i,j,k,nx,ny,nbl,nbr

        nx  = size(x)
        ny  = model % ny
        nbl = model % nbl
        nbr = model % nbr

        M = 0

        call model % coeff_L( x(1),  Y(iy(1,1):iy(1,ny)),   BL, MBL, ny, nbl )
        forall (i = 1:nbl)              A(ibl(i)) = BL(i)
        forall (i = 1:nbl, j = 1:ny)    M(ibl(i),iy(1,j)) = MBL(i,j)

        call model % coeff_R( x(nx), Y(iy(nx,1):iy(nx,ny)), BR, MBR, ny, nbr )
        forall (i = 1:nbr)              A(ibr(i)) = BR(i)
        forall (i = 1:nbr, j = 1:ny)    M(ibr(i),iy(nx,j)) = MBR(i,j)

        through_space: do i = 1, nx-1
            dx = x(i+1) - x(i)
            xm = (x(i+1) + x(i)) / 2

            forall (j = 1:ny)
                Ym(j) = ( Y(iy(i+1,j)) + Y(iy(i,j)) ) / 2
                Dm(j) = ( Y(iy(i+1,j)) - Y(iy(i,j)) ) / dx
            end forall

            call model % coeff(xm, Ym, Dm, Am, MY, MD, ny)

            ! Macierz( nr_rownynia, nr_niewiadomej )
            forall (j = 1:ny)   A(ia(i,j)) = Am(j)
            forall (j = 1:ny, k = 1:ny)
                M( ia(i,j), iy(i,  k) ) = MY(j,k) / 2 - MD(j,k) / dx
                M( ia(i,j), iy(i+1,k) ) = MY(j,k) / 2 + MD(j,k) / dx
            end forall
        end do through_space

    contains

        elemental integer function iy(ix,i)
            integer, intent(in) :: ix,i
            iy = ny*(ix-1) + i
        end function
        elemental integer function ia(ix,i)
            integer, intent(in) :: ix,i
            ia = nbl + i + ny*(ix-1)
        end function
        elemental integer function ibl(i)
            integer, intent(in) :: i
            ibl = i
        end function
        elemental integer function ibr(i)
            integer, intent(in) :: i
            ibr = nbl + ny*(nx-1) + i
        end function

    end subroutine

    subroutine model_advance(model, x, Y, M, dY)

        class(model_t) :: model
        real(fp), dimension(:), intent(in) :: x
        real(fp), dimension( size(x) * (model % ny) ), intent(in) :: Y
        real(fp), dimension( size(Y) ), intent(out) :: dY
        real(fp), dimension( size(dY), size(Y) ), intent(out)  :: M
        integer, dimension( size(dY) ) :: ipiv
        integer :: errno

        call model % matrix(x,Y,M,dY)

        dY = -dY
        call DGESV(size(M,2), 1, M, size(M,1), ipiv, dY, size(dY), errno)

        if ( errno > 0 ) then
            write (error_unit, '("Parameter ", I0, " had illegal value")') abs(errno)
            error stop "zero on matrix diagonyl: singular matrix"
        end if
        if ( errno < 0 ) then
            write (error_unit, '("Parameter ", I0, " had illegal value")') abs(errno)
            error stop "illegal parameter value"
        end if

    end subroutine

end module relaxation
