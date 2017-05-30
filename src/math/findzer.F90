module slf_findzer

    use iso_fortran_env
    use kind

    implicit none

    interface findzer
        module procedure findzer_single, findzer_array
    end interface


#   define FINDZER_ITER 128

contains

    subroutine findzer_single(x, xlo, xhi, delx, f)

        real(fp), intent(inout) :: x
        real(fp), intent(in) :: xlo,xhi,delx

        interface
            pure subroutine f(x,y,dy)
                import fp
                real(fp), intent(in) :: x
                real(fp), intent(out) :: y
                real(fp), intent(out), optional :: dy
            end subroutine
        end interface

        integer :: i
        real(fp) :: dx,yx,y,x1

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
        real(fp), intent(inout), dimension(nx) :: x
        real(fp), intent(in), dimension(nx)  :: xlo,xhi
        real(fp), intent(in)  :: delx

        interface
            pure subroutine f(x,y,dy)
                import fp
                real(fp), intent(in), dimension(:)  :: x
                real(fp), intent(out), dimension(size(x))  :: y
                real(fp), intent(out), dimension(size(x)) , optional :: dy
            end subroutine
        end interface

        integer :: i
        real(fp), dimension(nx) :: dx,yx,y,x1
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
