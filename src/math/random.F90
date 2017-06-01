module slf_random

    use iso_fortran_env

    implicit none

    character(len=*), parameter, private :: debugfmt = '(A,":",I0,": ",A)'
    real(real64), parameter, private :: pi = 4*atan(real(1,real64))

    interface random_number_gauss
        module subroutine random_number_gauss_one(a)
            real(real64), intent(out) :: a
        end subroutine

        module subroutine random_number_gauss_arr(A,n)
            real(real64), intent(out), dimension(n) :: A
            integer, intent(in) :: n
        end subroutine

        module subroutine random_number_gauss_2arr(A,B,n)
            real(real64), intent(out), dimension(n) :: A,B
            integer, intent(in) :: n
        end subroutine
    end interface random_number_gauss

contains

    subroutine random_init_urandom()
        integer :: n,err
        integer, allocatable :: seed(:)

        call random_seed(size = n)
        allocate(seed(n))
        open(33, file = '/dev/urandom', action = 'read', form = 'unformatted')
        read(33) seed
        close(33)
        call random_seed(put = seed)
        deallocate(seed)
    end subroutine

    subroutine random_number_gauss_one(a)
        real(real64), intent(out) :: a
        real(real64), dimension(2) :: U
        call random_number(U)
        a = sqrt( -2 * log(U(1)) ) * cos(2 * pi * U(2))
    end subroutine

    subroutine random_number_gauss_arr(A,n)
        real(real64), intent(out), dimension(n) :: A
        integer, intent(in) :: n
        real(real64), dimension(n,2) :: U
        call random_number(U)
        A = sqrt( -2 * log(U(:,1)) ) * cos(2 * pi * U(:,2))
    end subroutine

    subroutine random_number_gauss_2arr(A,B,n)
        real(real64), intent(out), dimension(n) :: A,B
        integer, intent(in) :: n
        real(real64), dimension(n,2) :: U
        real(real64), dimension(n) :: C
        call random_number(U)
        C = sqrt( -2 * log(U(:,1)) )
        A = C * cos(2 * pi * U(:,2))
        B = C * sin(2 * pi * U(:,2))
    end subroutine

end module
