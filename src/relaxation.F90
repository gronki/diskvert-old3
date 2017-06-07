module relaxation

    use iso_fortran_env

    use precision
    use slf_cgs
    use globals

    ! wspolczynniki wygenerowane
    use relax_coefficients

    implicit none

    interface
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
        procedure(coeff_driver_t), pointer, nopass :: get_AM => NULL()
        procedure(bound_driver_t), pointer, nopass :: get_BL => NULL()
        procedure(bound_driver_t), pointer, nopass :: get_BR => NULL()
        procedure(size_driver_t),  pointer, nopass :: get_sz => NULL()
    contains
        procedure :: matrix => model_matrix
    end type

    type(model_t), dimension(16), target, private :: model_collection


contains

    subroutine initialize_model_collection

        model_collection(1) % get_AM => COEFF_SS73DYF
        model_collection(1) % get_sz => COEFF_SS73DYF_SIZE
        model_collection(1) % get_BL => COEFF_SS73DYF_BL
        model_collection(1) % get_BR => COEFF_SS73DYF_BR
        model_collection(3) % get_AM => COEFF_MAGNDYF
        model_collection(3) % get_sz => COEFF_MAGNDYF_SIZE
        model_collection(3) % get_BL => COEFF_MAGNDYF_BL
        model_collection(3) % get_BR => COEFF_MAGNDYF_BR
        model_collection(2) % get_AM => COEFF_SS73COR
        model_collection(2) % get_sz => COEFF_SS73COR_SIZE
        model_collection(2) % get_BL => COEFF_SS73COR_BL
        model_collection(2) % get_BR => COEFF_SS73COR_BR
        model_collection(4) % get_AM => COEFF_MAGNCOR
        model_collection(4) % get_sz => COEFF_MAGNCOR_SIZE
        model_collection(4) % get_BL => COEFF_MAGNCOR_BL
        model_collection(4) % get_BR => COEFF_MAGNCOR_BR
        model_collection(8) % get_AM => COEFF_MAGNCORCND
        model_collection(8) % get_sz => COEFF_MAGNCORCND_SIZE
        model_collection(8) % get_BL => COEFF_MAGNCORCND_BL
        model_collection(8) % get_BR => COEFF_MAGNCORCND_BR

    end subroutine

    function model_nr(compton,magnetic,conduction)
        logical, intent(in) :: compton,magnetic,conduction
        integer :: model_nr

        model_nr = 1 + merge(1, 0, compton) + merge(2, 0, magnetic) + merge(4, 0, conduction)
    end function

    subroutine generate_coefficients_c_interface(X,nx,Y,ny,A,na,M) bind(C, name='generate_coefficients')

        integer(c_int), intent(in), value :: nx,ny,na
        real(c_double), intent(in), dimension(nx) :: X
        real(c_double), intent(in), dimension(nx*ny) :: Y
        real(c_double), intent(out), dimension(nx*na) :: A
        real(c_double), intent(out), dimension(nx*na,nx*ny) :: M

        type(model_t), pointer :: model

        call initialize_model_collection

        model => model_collection(model_nr( &
            & compton = .false., magnetic = .false., conduction = .false.))

        call model % matrix(x,Y,A,M)

    end subroutine

    subroutine model_matrix(model,x,Y,A,M)

        class(model_t) :: model
        integer :: ny, na, nbl, nbr
        real(fp), dimension(:), intent(in) :: x
        ! uklad: [ Y1(1) Y2(1) Y3(1) Y1(2) Y2(2) Y3(2) ... ]
        real(fp), dimension(:), target, intent(in) :: Y
        ! uklad: [ A1(1) A2(1) A3(1) A1(2) A2(2) A3(2) ... ]
        real(fp), dimension(:), target, intent(out) :: A
        real(fp), dimension(size(A),size(Y)), target, intent(out) :: M
        real(fp), dimension(:), allocatable :: Ym, Dm, Am
        real(fp), dimension(:,:), allocatable :: MY, MD
        real(fp), dimension(:), pointer :: BL, BR
        real(fp), dimension(:,:), pointer :: MBL, MBR
        real(fp) :: dx, xm
        integer :: i,nx

        nx = size(X)
        call model % get_sz(ny, na, nbl, nbr)

        if ( size(Y) /= nx*ny ) &
            & error stop "size of Y should be nx*ny"
        if ( size(A) /= nx*na - na + nbl + nbr ) &
            & error stop "size of Y should be (nx*na - na + nbl + nbr)"

        allocate( Ym(ny), Dm(ny), Am(na) )
        allocate( MY(na,ny), MD(na,ny) )

        BL  => A(1:nbl)
        MBL => M(1:nbl, 1:ny)
        BR  => A(nx*na - na + nbl : nx*na - na + nbl + nbr)
        MBR => M(nx*na - na + nbl : nx*na - na + nbl + nbr, 1 + ny*(nx-1) : ny*nx)

        call model % get_BL( x(1),  Y1(1),  BL, MBL, ny, nbl )
        call model % get_BR( x(nx), Y1(nx), BR, MBR, ny, nbr )

        through_space: do i = 1, nx-1
            dx = x(i+1) - x(i)
            xm = (x(i+1) + x(i)) / 2

            Ym = ( Y1(i+1) + Y1(i) ) / 2
            Dm = ( Y1(i+1) - Y1(i) ) / dx

            call model % get_AM(xm, Ym, Dm, Am, MY, MD, ny, na)

            ! Macierz( nr_rownania, nr_niewiadomej )
            A1(i) = Am
            M1(i,i)   = MY / 2 - MD / dx
            M1(i,i+1) = MY / 2 + MD / dx
        end do through_space

        deallocate( MY, MD, Ym, Dm, Am )

    contains

        function Y1(i)
            integer, intent(in) :: i
            real(fp), dimension(:), pointer :: Y1
            Y1 => Y(1 + ny*(i-1) : ny*i)
        end function

        function A1(i)
            integer, intent(in) :: i
            real(fp), dimension(:), pointer :: A1
            A1 => A(nbl + 1 + na*(i-1) : nbl + na*i)
        end function

        function M1(i,j)
            integer, intent(in) :: i,j
            real(fp), dimension(:,:), pointer :: M1
            M1 => M(nbl + 1 + na*(i-1) : nbl + na*i, 1 + ny*(j-1) : ny*j)
        end function


    end subroutine


end module relaxation
