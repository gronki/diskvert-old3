module slf_findzer

    use iso_fortran_env

    implicit none

    interface findzer
        module procedure findzer_single, findzer_array
    end interface


#   define FINDZER_ITER 128

contains

    subroutine findzer_single(x, xlo, xhi, delx, f)

        real(real64), intent(inout) :: x
        real(real64), intent(in) :: xlo,xhi,delx

        interface
            pure subroutine f(x,y,dy)
                import real64
                real(real64), intent(in) :: x
                real(real64), intent(out) :: y
                real(real64), intent(out), optional :: dy
            end subroutine
        end interface

        integer :: i
        real(real64) :: dx,yx,y,x1

        main_loop : do i = 1,FINDZER_ITER

            call f(x,y,yx)
            if (yx /= 0) then
                dx = - y / yx
            elseif ( dx == 0 ) then
                exit main_loop
            end if

            x1 = x + dx
            if ( x1 > xhi ) x1 = (x + xhi) / 2
            if ( x1 < xlo ) x1 = (x + xlo) / 2
            x = x1

            if (abs(dx) < abs(delx)) exit main_loop

        end do main_loop

    end subroutine

    subroutine findzer_array(x, nx, xlo, xhi, delx, f)

        integer, intent(in) :: nx
        real(real64), intent(inout), dimension(nx) :: x
        real(real64), intent(in), dimension(nx)  :: xlo,xhi
        real(real64), intent(in)  :: delx

        interface
            pure subroutine f(x,y,dy)
                import real64
                real(real64), intent(in), dimension(:)  :: x
                real(real64), intent(out), dimension(size(x))  :: y
                real(real64), intent(out), dimension(size(x)) , optional :: dy
            end subroutine
        end interface

        integer :: i
        real(real64), dimension(nx) :: dx,yx,y,x1
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

                x1 = x + dx
                where ( x1 > xhi ) x1 = (x + xhi) / 2
                where ( x1 < xlo ) x1 = (x + xlo) / 2
                x = x1

                converged = abs(dx) < abs(delx)
            end where

            if (all(converged)) exit main_loop

        end do main_loop

    end subroutine

end module
