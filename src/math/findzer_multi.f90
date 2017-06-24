module slf_findzermulti

    use iso_fortran_env
    use ieee_arithmetic
    
    use precision

    implicit none

contains

    subroutine findzer_multi(x, nx, xlo, xhi, delx, nb, f)

        real(fp), intent(inout) :: x(:)
        real(fp), intent(in) :: xlo,xhi,delx
        integer, intent(in) :: nb
        integer, intent(out) :: nx

        interface
            pure subroutine f(x,y,dy)
                import fp
                real(fp), intent(in) :: x
                real(fp), intent(out) :: y
                real(fp), intent(out), optional :: dy
            end subroutine
        end interface

        real(fp) :: y(nb), dy(nb), x0(nb)
        integer :: i


        do concurrent (i=1:nb)
            x0(i) = (i-1)/(nb-1.) * (xhi - xlo) + xlo
        end do

        call f(x0(1),y(1))
        nx = 0
        scan_for_solutions: do i = 2, nb
            call f(x0(i),y(i))
            if ( y(i-1) * y(i) < 0 ) then
                call linsect(x(nx+1), x0(i-1), x0(i), delx)
                nx = nx + 1
            end if
            if ( nx >= size(x) ) exit scan_for_solutions
        end do scan_for_solutions

    contains

        subroutine linsect(x,xlo0,xhi0,delx)
            real(fp), intent(in) :: xhi0,xlo0,delx
            real(fp), intent(out) :: x
            real(fp) :: xhi,xlo,yhi,ylo,y,s
            integer :: i,n

            n = log(abs((xhi0-xlo0)/delx))/log(2.0)

            xlo = xlo0
            xhi = xhi0

            sect_loop: do i = 1, n
                call f(xlo,ylo)
                call f(xhi,yhi)
                s = sign(real(1,kind=fp),yhi-ylo)

                x = ( yhi * xlo - ylo * xhi ) / ( yhi - ylo )
                call f(x,y)

                if ( y * s > 0 ) then
                    xhi = x
                elseif ( y * s < 0 ) then
                    xlo = x
                else
                    exit sect_loop
                end if

                if ( abs(xhi-xlo) < abs(delx) ) exit sect_loop
            end do sect_loop

        end subroutine

    end subroutine

end module
