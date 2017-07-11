module slf_random

    use iso_fortran_env, only: r64 => real64
    implicit none

    character(len=*), parameter, private :: debugfmt = '(A,":",I0,": ",A)'
    real(r64), parameter, private :: pi = 4*atan(real(1,r64))

    interface random_number_gauss
#     if __GNUC__ < 6
        module procedure :: random_number_gauss_one
        module procedure :: random_number_gauss_arr
        module procedure :: random_number_gauss_2arr
#     else
        module subroutine random_number_gauss_one(a)
          real(r64), intent(out) :: a
        end subroutine

        module subroutine random_number_gauss_arr(B)
          real(r64), intent(out), dimension(:) :: B
        end subroutine

        module subroutine random_number_gauss_2arr(C,D)
          real(r64), intent(out), dimension(:) :: C
          real(r64), intent(out), dimension(size(C)) :: D
        end subroutine
#     endif
    end interface random_number_gauss

contains

    subroutine random_seed_urandom()
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
        real(r64), intent(out) :: a
        real(r64), dimension(2) :: U
        call random_number(U)
        a = sqrt( -2 * log(U(1)) ) * cos(2 * pi * U(2))
    end subroutine

    subroutine random_number_gauss_arr(A)
        real(r64), intent(out), dimension(:) :: A
        real(r64), dimension(size(A),2) :: U
        call random_number(U)
        A = sqrt( -2 * log(U(:,1)) ) * cos(2 * pi * U(:,2))
    end subroutine

    subroutine random_number_gauss_2arr(A,B)
        real(r64), intent(out), dimension(:) :: A
        real(r64), intent(out), dimension(size(A)) :: B
        real(r64), dimension(size(A),2) :: U
        real(r64), dimension(size(A)) :: C
        call random_number(U)
        C = sqrt( -2 * log(U(:,1)) )
        A = C * cos(2 * pi * U(:,2))
        B = C * sin(2 * pi * U(:,2))
    end subroutine

end module
