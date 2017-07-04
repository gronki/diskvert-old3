module slf_findzer

    use iso_fortran_env, only: r64 => real64
    implicit none

    interface findzer
        module procedure findzer_single, findzer_array
    end interface

    integer, parameter :: findzer_iter = 64
    private :: ramp

contains

    subroutine findzer_single(x, xlo, xhi, delx, f)

        real(r64), intent(inout) :: x
        real(r64), intent(in) :: xlo,xhi,delx

        interface
            pure subroutine f(x,y,dy)
                import r64
                real(r64), intent(in) :: x
                real(r64), intent(out) :: y
                real(r64), intent(out), optional :: dy
            end subroutine
        end interface

        integer :: i
        real(r64) :: dx,yx,y,x1

        main_loop : do i = 1,FINDZER_ITER

            call f(x,y,yx)
            if (yx /= 0) then
                dx = - y / yx
            elseif ( dx == 0 ) then
                exit main_loop
            end if

            x1 = x + dx * ramp(i,findzer_iter)
            if ( x1 > xhi ) x1 = (x + xhi) / 2
            if ( x1 < xlo ) x1 = (x + xlo) / 2
            x = x1

            if (abs(dx) < abs(delx)) exit main_loop

        end do main_loop

    end subroutine

    subroutine findzer_array(x, nx, xlo, xhi, delx, f)

        integer, intent(in) :: nx
        real(r64), intent(inout), dimension(nx) :: x
        real(r64), intent(in), dimension(nx)  :: xlo,xhi
        real(r64), intent(in)  :: delx

        interface
            pure subroutine f(x,y,dy)
                import r64
                real(r64), intent(in), dimension(:)  :: x
                real(r64), intent(out), dimension(size(x))  :: y
                real(r64), intent(out), dimension(size(x)) , optional :: dy
            end subroutine
        end interface

        integer :: i
        real(r64), dimension(nx) :: dx,yx,y,x1
        logical, dimension(nx) :: converged

        converged = .false.
        dx = 0

        main_loop : do i = 1,FINDZER_ITER

            call f(x,y,yx)

            where ( .not. converged )

                where ( yx /= 0 )
                    dx = - y / yx
                elsewhere ( dx == 0 )
                    converged = .true.
                end where

                x1 = x + dx * ramp(i,findzer_iter)
                where ( x1 > xhi ) x1 = (x + xhi) / 2
                where ( x1 < xlo ) x1 = (x + xlo) / 2
                x = x1

                converged = abs(dx) < abs(delx)
            end where

            if (all(converged)) exit main_loop

        end do main_loop

    end subroutine

    elemental function ramp(i,n) result(y)
      integer, intent(in) :: i,n
      real(r64) :: t,y
      t = merge(real(i) / n, 1.0, i .lt. n)
      y = (3 - 2*t) * t**2
    end function

end module
