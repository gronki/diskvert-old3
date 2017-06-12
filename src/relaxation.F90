module relaxation

    use iso_fortran_env

    use precision
    use slf_cgs
    use globals
    use lapack

    ! wspolczynniki wygenerowane
    use relax_coefficients

    implicit none

    abstract interface
        pure subroutine coeff_driver_t(z, Y, D, A, MY, MD, ny, na)
            import fp
            integer, intent(in) :: ny,na
            real(fp), intent(in) :: z
            real(fp), intent(in), dimension(ny) :: Y
            real(fp), intent(in), dimension(ny) :: D
            real(fp), intent(out), dimension(na) :: A
            real(fp), intent(out), dimension(na,ny) :: MY
            real(fp), intent(out), dimension(na,ny) :: MD
        end subroutine
        pure subroutine bound_driver_t(z, Y, B, M, ny, nb)
            import fp
            integer, intent(in) :: ny,nb
            real(fp), intent(in) :: z
            real(fp), intent(in), dimension(ny) :: Y
            real(fp), intent(out), dimension(nb) :: B
            real(fp), intent(out), dimension(nb,ny) :: M
        end subroutine
        pure subroutine size_driver_t(ny, neq, nbl, nbr)
            integer, intent(out) :: ny
            integer, intent(out) :: neq
            integer, intent(out) :: nbl
            integer, intent(out) :: nbr
        end subroutine size_driver_t
    end interface

    type :: model_t
        logical, private :: initialized = .FALSE.
        procedure(coeff_driver_t), pointer, nopass, private :: get_AM
        procedure(bound_driver_t), pointer, nopass, private :: get_BL, get_BR
        procedure(size_driver_t),  pointer, nopass, private :: get_sz
        integer, private :: ny, na, nbl, nbr
        integer, private :: nl, nu
    contains
        procedure :: matrix  => model_matrix
        procedure :: init    => model_initialize
        procedure :: advance => model_corrections
    end type

contains

    subroutine model_initialize (model, compton, magnetic, conduction)
        class(model_t) :: model
        logical, intent(in) :: compton,magnetic,conduction
        integer :: model_nr

        model_nr = 1 + merge(1, 0, compton) + merge(2, 0, magnetic) &
            & + merge(4, 0, conduction)

        select case(model_nr)
        case(1)
            model % get_AM => COEFF_SS73DYF
            model % get_sz => COEFF_SS73DYF_SIZE
            model % get_BL => COEFF_SS73DYF_BL
            model % get_BR => COEFF_SS73DYF_BR
        case(2)
            model % get_AM => COEFF_SS73COR
            model % get_sz => COEFF_SS73COR_SIZE
            model % get_BL => COEFF_SS73COR_BL
            model % get_BR => COEFF_SS73COR_BR
        case(3)
            model % get_AM => COEFF_MAGNDYF
            model % get_sz => COEFF_MAGNDYF_SIZE
            model % get_BL => COEFF_MAGNDYF_BL
            model % get_BR => COEFF_MAGNDYF_BR
        case(4)
            model % get_AM => COEFF_MAGNCOR
            model % get_sz => COEFF_MAGNCOR_SIZE
            model % get_BL => COEFF_MAGNCOR_BL
            model % get_BR => COEFF_MAGNCOR_BR
        case(5)
            model % get_AM => COEFF_SS73DYFCND
            model % get_sz => COEFF_SS73DYFCND_SIZE
            model % get_BL => COEFF_SS73DYFCND_BL
            model % get_BR => COEFF_SS73DYFCND_BR
        case(6)
            model % get_AM => COEFF_SS73CORCND
            model % get_sz => COEFF_SS73CORCND_SIZE
            model % get_BL => COEFF_SS73CORCND_BL
            model % get_BR => COEFF_SS73CORCND_BR
        case(7)
            model % get_AM => COEFF_MAGNDYFCND
            model % get_sz => COEFF_MAGNDYFCND_SIZE
            model % get_BL => COEFF_MAGNDYFCND_BL
            model % get_BR => COEFF_MAGNDYFCND_BR
        case(8)
            model % get_AM => COEFF_MAGNCORCND
            model % get_sz => COEFF_MAGNCORCND_SIZE
            model % get_BL => COEFF_MAGNCORCND_BL
            model % get_BR => COEFF_MAGNCORCND_BR
        case default
            error stop "This model is not implemented."
        end select

        call model % get_sz(    model % ny,     &
                            &   model % na,     &
                            &   model % nbl,    &
                            &   model % nbr)

        model % nl = (model % na) + (model % nbl) - 1
        model % nu = 2 * (model % ny) - (model % nbl) - 1
        model % initialized = .TRUE.

    end subroutine


    subroutine model_matrix(model,x,Y,A,M)

        class(model_t) :: model
        real(fp), dimension(:), intent(in) :: x
        ! uklad: [ Y1(1) Y2(1) Y3(1) Y1(2) Y2(2) Y3(2) ... ]
        real(fp), dimension(size(x) * (model % ny)), intent(in) :: Y
        ! uklad: [ A1(1) A2(1) A3(1) A1(2) A2(2) A3(2) ... ]
        real(fp), dimension(size(x) * (model % na)), intent(out) :: A
        real(fp), dimension(size(A), size(Y)), target, intent(out) :: M

        real(fp), dimension(model % ny) :: Ym, Dm
        real(fp), dimension(model % na) :: Am
        real(fp), dimension(model % na, model % ny) :: MY, MD
        real(fp), dimension(model % nbl) :: BL
        real(fp), dimension(model % nbl, model % ny) :: MBL
        real(fp), dimension(model % nbr) :: BR
        real(fp), dimension(model % nbr, model % ny) :: MBR

        real(fp) :: dx, xm
        integer :: i,j,k,nx,ny,na,nbl,nbr

        if ( .NOT. model % initialized ) then
            error stop "model not initialized"
        end if


        nx  = size(X)
        ny  = model % ny
        na  = model % na
        nbl = model % nbl
        nbr = model % nbr

        M = 0

        call model % get_BL( x(1),  Y(iy(1,1):iy(1,ny)),   BL, MBL, ny, nbl )
        forall (i = 1:nbl)              A(ibl(i)) = BL(i)
        forall (i = 1:nbl, j = 1:ny)    M(ibl(i),iy(1,j)) = MBL(i,j)

        call model % get_BR( x(nx), Y(iy(nx,1):iy(nx,ny)), BR, MBR, ny, nbr )
        forall (i = 1:nbr)              A(ibr(i)) = BR(i)
        forall (i = 1:nbr, j = 1:ny)    M(ibr(i),iy(nx,j)) = MBR(i,j)

        through_space: do i = 1, nx-1
            dx = x(i+1) - x(i)
            xm = (x(i+1) + x(i)) / 2

            forall (j = 1:ny)
                Ym(j) = ( Y(iy(i+1,j)) + Y(iy(i,j)) ) / 2
                Dm(j) = ( Y(iy(i+1,j)) - Y(iy(i,j)) ) / dx
            end forall

            call model % get_AM(xm, Ym, Dm, Am, MY, MD, ny, na)

            ! Macierz( nr_rownania, nr_niewiadomej )
            forall (j = 1:na)   A(ia(i,j)) = Am(j)
            forall (j = 1:na, k = 1:ny)
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
            ia = nbl + i + na*(ix-1)
        end function
        elemental integer function ibl(i)
            integer, intent(in) :: i
            ibl = i
        end function
        elemental integer function ibr(i)
            integer, intent(in) :: i
            ibr = nbl + na*(nx-1) + i
        end function


    end subroutine

    subroutine model_corrections(model, x, Y, A, M, dY)

        class(model_t) :: model
        real(fp), dimension(:), intent(in) :: x
        real(fp), dimension( size(x) * (model % ny) ), intent(in) :: Y
        real(fp), dimension( size(x) * (model % na) ), intent(out) :: A
        real(fp), dimension( size(A), size(Y) ), intent(out)  :: M
        real(fp), dimension( size(Y) ), intent(out) :: dY
        integer, dimension( size(A) ) :: ipiv
        integer :: errno

        if ( .not. model % initialized ) error stop "model not initialized"

        call model % matrix(x,Y,A,M)

        dY = - A
        call DGESV(size(M,2), 1, M, size(M,1), ipiv, dY, size(A), errno)

        if ( errno > 0 ) then
            write (error_unit, '("Parameter ", I0, " had illegal value")') abs(errno)
            error stop "zero on matrix diagonal: singular matrix"
        end if
        if ( errno < 0 ) then
            write (error_unit, '("Parameter ", I0, " had illegal value")') abs(errno)
            error stop "illegal parameter value"
        end if

    end subroutine



end module relaxation
