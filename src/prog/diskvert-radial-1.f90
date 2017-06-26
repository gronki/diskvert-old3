program diskvert_radial_1

    use precision
    use slf_cgs
    use diskapprox
    use globals
    use settings
    use setup

    implicit none

    integer :: errno,i,j
    integer, parameter :: n = 2**8
    real(fp), parameter :: r0 = 3.1
    real(fp) :: r12, r23
    real(fp), dimension(n) :: r, rho, T, H, rho2, T2, H2
    real(fp) :: M(3,3), A(3)
    integer :: ipiv(3)

    errno = 0
    call init(errno, confread)

    call apx_zonebounds(r12,r23)

    forall (i = 1:n)
        r(i) = r0 * exp( (i-1)/real(n-1) * (log(r23*10)-log(r0)) )
    end forall

    call apx_sel(r, r12, r23, rho, T, H)

    rho2 = rho
    T2 = T
    H2 = H

    compute_precise_model: do i = 1,n
        converge: do j = 1,16
            call apx_matrix(r(i), rho2(i), T2(i), H2(i), A, M)
            call dgesv(3, 1, M, 3, ipiv, A, 3, errno)
            rho2(i) = rho2(i) + A(1) * ramp(j,16)
            T2(i) = T2(i) + A(2) * ramp(j,16)
            H2(i) = H2(i) + A(3) * ramp(j,16)
        end do converge
    end do compute_precise_model

    open(33, file = 'profile.dat', action = 'write')

    do i = 1,n
        write (33, '(I7,7Es12.4)') i, r(i), rho(i), rho2(i), T(i), T2(i), &
                H(i), H2(i)
    end do

    close(33)

contains

    elemental function ramp(i,n) result(y)
        integer, intent(in) :: i,n
        real(fp) :: t,y
        t = merge(real(i) / n, 1.0, i .le. n)
        y = 3 * t**2 - 2 * t**3
    end function

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
