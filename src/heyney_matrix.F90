module heyney_matrix

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
        ! model_collection(1) % get_BL => COEFF_SS73DYF_BL
        ! model_collection(1) % get_BR => COEFF_SS73DYF_BR
        model_collection(3) % get_AM => COEFF_MAGNDYF
        model_collection(3) % get_sz => COEFF_MAGNDYF_SIZE
        ! model_collection(3) % get_BL => COEFF_MAGNDYF_BL
        ! model_collection(3) % get_BR => COEFF_MAGNDYF_BR
        model_collection(2) % get_AM => COEFF_SS73COR
        model_collection(2) % get_sz => COEFF_SS73COR_SIZE
        ! model_collection(2) % get_BL => COEFF_SS73COR_BL
        ! model_collection(2) % get_BR => COEFF_SS73COR_BR
        model_collection(4) % get_AM => COEFF_MAGNCOR
        model_collection(4) % get_sz => COEFF_MAGNCOR_SIZE
        ! model_collection(4) % get_BL => COEFF_MAGNCOR_BL
        ! model_collection(4) % get_BR => COEFF_MAGNCOR_BR
        model_collection(8) % get_AM => COEFF_MAGNCORCND
        model_collection(8) % get_sz => COEFF_MAGNCORCND_SIZE
        ! model_collection(8) % get_BL => COEFF_MAGNCORCND_BL
        ! model_collection(8) % get_BR => COEFF_MAGNCORCND_BR

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
        real(fp), dimension(:), intent(in) :: Y
        ! uklad: [ A1(1) A2(1) A3(1) A1(2) A2(2) A3(2) ... ]
        real(fp), dimension(:), intent(out) :: A
        real(fp), dimension(size(A),size(Y)), intent(out) :: M
        real(fp), dimension(:), allocatable :: Am, Ym, Dm
        real(fp), dimension(:,:), allocatable :: MY, MD
        real(fp) :: dx, xm
        integer :: i

        call model % get_sz(ny, na, nbl, nbr)

        allocate( Am(na), Ym(ny), Dm(ny) )
        allocate( MY(na,ny), MD(na,ny) )

        through_space: do i = 1, size(x)-1
            dx = x(i+1) - x(i)
            xm = (x(i+1) + x(i)) / 2

            Ym = (Y(ilo(i):ihi(i)) + Y(ilo(i+1):ihi(i+1))) / 2
            Dm = (Y(ilo(i+1):ihi(i+1)) - Y(ilo(i):ihi(i))) / dx

            call model % get_AM(xm, Ym, Dm, Am, MY, MD, ny, na)

            A(ilo(i) : ihi(i)) = Am
            ! Macierz( nr_rownania, nr_niewiadomej )
            M(ilo(i) : ihi(i), ilo(i) : ihi(i)) = MY / 2 - MD / dx
            M(ilo(i) : ihi(i), ilo(i+1) : ihi(i+1)) = MY / 2 + MD / dx
        end do through_space

        deallocate( Am, Ym, Dm, MY, MD )

    contains

        elemental integer function ilo(i)
            integer, intent(in) :: i
            ilo = 1 + ny*(i-1)
        end function
        elemental integer function ihi(i)
            integer, intent(in) :: i
            ihi = na*i
        end function

    end subroutine


end module heyney_matrix
