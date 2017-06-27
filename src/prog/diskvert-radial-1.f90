program diskvert_radial_1

    use precision
    use slf_cgs
    use diskapprox
    use globals
    use settings
    use setup
    use model_ss73

    implicit none

    integer :: errno,i,j
    integer, parameter :: n = 2**9
    real(fp) :: r12, r23
    real(fp), dimension(n) :: r
    real(fp), dimension(n) :: rho1, T1, H1
    real(fp), dimension(n) :: rho2, T2, H2
    real(fp), dimension(n) :: rho3, T3, H3
    real(fp) :: M(3,3), A(3)
    integer :: ipiv(3)

    integer, parameter :: nz = 2**11
    real(fp) :: z(nz), y(ny,nz), dy(ny,nz), ap(na,nz)
    integer :: nmax

    errno = 0
    call init(errno, confread)

    call apx_zonebounds(r12,r23)

    forall (i = 1:n)
        r(i) = 3 * exp( i / real(n) * (log(r23*10)-log(3d0)) )
    end forall

    call apx_sel(r, r12, r23, rho1, T1, H1)

    rho2 = rho1
    T2 = T1
    H2 = H1

    compute_precise_model: do i = 1,n
        ! converge: do j = 1,16
        !     call apx_matrix(r(i), rho2(i), T2(i), H2(i), A, M)
        !     call dgesv(3, 1, M, 3, ipiv, A, 3, errno)
        !     rho2(i) = rho2(i) + A(1) * ramp(j,16)
        !     T2(i) = T2(i) + A(2) * ramp(j,16)
        !     H2(i) = H2(i) + A(3) * ramp(j,16)
        ! end do converge
        call diskapx2(r(i), rho2(i), T2(i), H2(i))
        r_calc = r(i)
        call eval_globals
        call run_ss73(z,nz,y,dy,ap,nmax)
        rho3(i) = ap(c_rho,nz)
        T3(i) = ap(c_Tgas,nz)
        H3(i) = z(1) / 2
        search_height: do j = 1,nz
            if (ap(c_rho,j) .ge. rho3(i) / 2) then
                H3(i) = z(j)
                exit search_height
            end if
        end do search_height
    end do compute_precise_model

    do i = 1,n
        write (6, '(I7,10(X,Es11.4))') i, r(i), rho1(i), T1(i), H1(i), &
            rho2(i), T2(i), H2(i), rho3(i), T3(i), H3(i)
    end do

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
