module slf_deriv

    use iso_fortran_env, only: r64 => real64
    implicit none

contains

    subroutine deriv(x,yin,yout)
        real(r64), intent(in) :: x(:), yin(size(x))
        real(r64), intent(out) :: yout(size(x))
        real(r64) :: dx2
        integer :: i,n

        n = size(yout)

        if ( n .lt. 3 ) then
            write (0,*) 'n should be at least 3!'
            return
        end if

        do i=2,n-1
            yout(i) = (1d0*yin(i+1)-yin(i-1))/(1d0*x(i+1)-x(i-1))
        end do
        dx2 = 2 * ( yin(3) + yin(1) - 2 * yin(2) ) / ( (x(3)-x(2))**2 + (x(2)-x(1))**2 )
        yout(1) = yout(2) + (x(1)-x(2))*dx2
        dx2 = 2 * ( yin(n) + yin(n-2) - 2 * yin(n-1) ) / ( (x(n)-x(n-1))**2 + (x(n-1)-x(n-2))**2 )
        yout(n) = yout(n-1) + (x(n)-x(n-1))*dx2
    end subroutine

    subroutine deriv2(x,yin,yout)
        real(r64), intent(in) :: x(:), yin(size(x))
        real(r64), intent(out) :: yout(size(x))
        integer :: i,n
        n = size(yout)

        do i=2,n-1
            yout(i) = 2 * ( yin(i+1) + yin(i-1) - 2 * yin(i) ) / ( (x(i+1)-x(i))**2 + (x(i-1)-x(i))**2 )
        end do
        yout(1) = yout(2) + ( x(1) - x(2) )*( yout(3)-yout(2) )/( x(3)-x(2) )
        yout(n) = yout(n-1) + ( x(n) - x(n-1) )*( yout(n-1)-yout(n-2) )/( x(n-1)-x(n-2) )
    end subroutine


end module
