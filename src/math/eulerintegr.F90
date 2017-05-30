module slf_eulerintegr

    use ieee_arithmetic
    use iso_fortran_env

    implicit none

    interface eulerintegr
        module procedure eulerintegr_full
        module procedure eulerintegr_simple
    end interface eulerintegr

contains

    subroutine eulerintegr_full(x, y0, fder, na, y, dy, a, nmax)

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

        real(real64) :: dx
        logical :: abort
        integer :: i

        y(:,1) = y0
        a(:,1) = 0
        dy(:,1) = 0

        abort = .false.

        integr_loop : do i = 1,size(x)

            if ( i > 1 ) then
                a(:,i) = a(:,i-1)
                dy(:,i) = dy(:,i-1)
            end if

            call fder(x(i), y(:,i), dy(:,i), a(:,i), abort)

            if ( i < size(x) ) then
                dx = x(i+1) - x(i)
                y(:,i+1) = y(:,i) + dx * dy(:,i)
            end if

            if ( abort ) exit integr_loop
            if ( any(.not.ieee_is_normal(dy(:,i))) ) exit integr_loop

            if ( present(nmax) ) nmax = i
        end do integr_loop

    end subroutine

    subroutine eulerintegr_simple(x, y0, fder, y, dy, nmax)

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

        real(real64) :: dx
        logical :: abort
        integer :: i

        y(:,1) = y0
        dy(:,1) = 0

        abort = .false.

        integr_loop : do i = 1,size(x)

            if ( i > 1 ) then
                dy(:,i) = dy(:,i-1)
            end if

            call fder(x(i), y(:,i), dy(:,i), abort)

            if ( i < size(x) ) then
                dx = x(i+1) - x(i)
                y(:,i+1) = y(:,i) + dx * dy(:,i)
            end if

            if ( abort ) exit integr_loop
            if ( any(.not.ieee_is_normal(dy(:,i))) ) exit integr_loop

            if ( present(nmax) ) nmax = i
        end do integr_loop

    end subroutine

end module
