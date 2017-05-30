module slf_rk4integr

    use iso_fortran_env
    use ieee_arithmetic

    implicit none

    interface rk4integr
        module procedure rk4integr_full
        module procedure rk4integr_simple
    end interface

contains

    subroutine rk4integr_full(x, y0, fder, na, y, dy, a, nmax)

        real(real64), intent(in) :: x(:), y0(:)
        integer, intent(in) :: na
        real(real64), intent(out) :: y(size(y0),size(x))
        real(real64), intent(out) :: dy(size(y0),size(x))
        real(real64), intent(out) :: a(na,size(x))
        integer, intent(out), optional :: nmax

        interface
            subroutine fder(x,y,dy,a,abort)
                import real64
                real(real64), intent(in) :: x, y(:)
                real(real64), intent(inout) :: dy(size(y)), a(:)
                logical, intent(inout) :: abort
            end subroutine
        end interface

        real(real64) :: dx, dy1(size(y0)), dy2(size(y0)), dy3(size(y0)), dy4(size(y0))
        logical :: abort
        integer :: i

        y(:,1) = y0
        a(:,1) = 0
        dy(:,1) = 0

        abort = .false.
        nmax = 0

        integr_loop : do i = 1,size(x)

            if ( i > 1 ) then
                a(:,i) = a(:,i-1)
                dy(:,i) = dy(:,i-1)
            end if

            call fder(x(i), y(:,i), dy1, a(:,i), abort)
            if (abort) exit integr_loop

            if ( i < size(x) ) then
                dx = x(i+1) - x(i)
                a(:,i+1) = a(:,i)

                call fder(x(i) + (dx/2), y(:,i) + dy1 * (dx/2), dy2, a(:,i+1), abort)
                if (abort) exit integr_loop
                call fder(x(i) + (dx/2), y(:,i) + dy2 * (dx/2), dy3, a(:,i+1), abort)
                if (abort) exit integr_loop
                call fder(x(i+1), y(:,i) + dy3 * dx, dy4, a(:,i+1), abort)
                if (abort) exit integr_loop

                dy(:,i) = ( dy1 + 2*dy2 + 2*dy3 + dy4 ) / 6
                if ( any(.not.ieee_is_normal(dy(:,i))) ) exit integr_loop

                y(:,i+1) = y(:,i) + dx * dy(:,i)
            else
                dy(:,i) = dy1
            end if

            if (present(nmax)) nmax = i

        end do integr_loop

    end subroutine

    subroutine rk4integr_simple(x, y0, fder, y, dy, nmax)

        real(real64), intent(in) :: x(:), y0(:)
        real(real64), intent(out) :: y(size(y0),size(x))
        real(real64), intent(out) :: dy(size(y0),size(x))
        integer, intent(out), optional :: nmax

        interface
            subroutine fder(x,y,dy,abort)
                import real64
                real(real64), intent(in) :: x, y(:)
                real(real64), intent(inout) :: dy(size(y))
                logical, intent(inout) :: abort
            end subroutine
        end interface

        real(real64) :: dx, dy1(size(y0)), dy2(size(y0)), dy3(size(y0)), dy4(size(y0))
        logical :: abort
        integer :: i

        y(:,1) = y0
        dy(:,1) = 0

        abort = .false.
        nmax = 0

        integr_loop : do i = 1,size(x)

            if ( i > 1 ) then
                dy(:,i) = dy(:,i-1)
            end if

            call fder(x(i), y(:,i), dy1, abort)
            if (abort) exit integr_loop

            if ( i < size(x) ) then
                dx = x(i+1) - x(i)

                call fder(x(i) + (dx/2), y(:,i) + dy1 * (dx/2), dy2, abort)
                if (abort) exit integr_loop
                call fder(x(i) + (dx/2), y(:,i) + dy2 * (dx/2), dy3, abort)
                if (abort) exit integr_loop
                call fder(x(i+1), y(:,i) + dy3 * dx, dy4, abort)
                if (abort) exit integr_loop

                dy(:,i) = ( dy1 + 2*dy2 + 2*dy3 + dy4 ) / 6
                if ( any(.not.ieee_is_normal(dy(:,i))) ) exit integr_loop

                y(:,i+1) = y(:,i) + dx * dy(:,i)
            else
                dy(:,i) = dy1
            end if

            if (present(nmax)) nmax = i

        end do integr_loop

    end subroutine

end module
