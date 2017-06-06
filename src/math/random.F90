module slf_random

    use iso_fortran_env
    
    use precision

    implicit none

    character(len=*), parameter, private :: debugfmt = '(A,":",I0,": ",A)'
    real(fp), parameter, private :: pi = 4*atan(real(1,fp))

    interface random_number_gauss
        module subroutine random_number_gauss_one(a)
            real(fp), intent(out) :: a
        end subroutine

        module subroutine random_number_gauss_arr(B)
            real(fp), intent(out), dimension(:) :: B
        end subroutine

        module subroutine random_number_gauss_2arr(C,D)
            real(fp), intent(out), dimension(:) :: C
            real(fp), intent(out), dimension(size(C)) :: D
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
        real(fp), intent(out) :: a
        real(fp), dimension(2) :: U
        call random_number(U)
        a = sqrt( -2 * log(U(1)) ) * cos(2 * pi * U(2))
    end subroutine

    subroutine random_number_gauss_arr(A)
        real(fp), intent(out), dimension(:) :: A
        real(fp), dimension(size(A),2) :: U
        call random_number(U)
        A = sqrt( -2 * log(U(:,1)) ) * cos(2 * pi * U(:,2))
    end subroutine

    subroutine random_number_gauss_2arr(A,B)
        real(fp), intent(out), dimension(:) :: A
        real(fp), intent(out), dimension(size(A)) :: B
        real(fp), dimension(size(A),2) :: U
        real(fp), dimension(size(A)) :: C
        call random_number(U)
        C = sqrt( -2 * log(U(:,1)) )
        A = C * cos(2 * pi * U(:,2))
        B = C * sin(2 * pi * U(:,2))
    end subroutine

end module
