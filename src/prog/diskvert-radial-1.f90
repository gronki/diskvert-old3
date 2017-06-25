program diskvert_radial_1

    use precision
    use slf_cgs
    use diskapprox
    use globals
    use settings
    use setup

    implicit none

    integer :: errno,i
    integer, parameter :: n = 2**10
    real(fp), parameter :: r0 = 3.1
    real(fp) :: r12, r23
    real(fp), dimension(n) :: r, rho, T, H

    errno = 0
    call init(errno, confread)

    call apx_zonebounds(r12,r23,24)

    forall (i = 1:n)
        r(i) = r0 * exp( (i-1)/real(n-1) * (log(r23*5)-log(r0)) )
    end forall

    call apx_sel(r, r12, r23, rho, T, H)

    open(33, file = 'profile.dat', action = 'write')

    do i = 1,n
        write (33, '(I7,4Es12.4)') i, r(i), rho(i), T(i), H(i)
    end do

    close(33)

contains

    subroutine confread(cfg,errno)
        type(config), intent(inout) :: cfg
        integer, intent(inout) :: errno
        character(128) :: buf

        if ( .not. cfg % contains('alpha') ) then
            error stop "alpha key is required!"
        end if

        call cfg%get('alpha',buf)
        read (buf,*) alpha
    end subroutine

end program
