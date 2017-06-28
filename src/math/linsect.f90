module slf_linsect

    use iso_fortran_env, only: r64 => real64
    implicit none

contains

    subroutine linsect(f,x,xlo0,xhi0,delx)
        real(r64), intent(in) :: xhi0,xlo0,delx
        real(r64), intent(out) :: x
        real(r64) :: xhi,xlo,yhi,ylo,y,s
        integer :: i,n

        interface
            pure subroutine f(x,y)
                import r64
                real(r64), intent(in) :: x
                real(r64), intent(out) :: y
            end subroutine
        end interface

        n = log(abs((xhi0-xlo0)/delx))/log(2.0)

        xlo = xlo0
        xhi = xhi0

        sect_loop: do i = 1, n
            call f(xlo,ylo)
            call f(xhi,yhi)
            s = sign(real(1,r64),yhi-ylo)

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

end module
