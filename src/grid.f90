module grid

    use iso_fortran_env, only: r64 => real64
    implicit none

contains

    elemental function space_linear(i,n,h) result(f)
        integer, intent(in) :: i,n
        real(r64), intent(in) :: h
        real(r64) :: f,x
        x = real(i - 1, r64) / real(n - 1, r64)
        f = x * h
    end function

    elemental function space_linear_r(i,n,h) result(f)
        integer, intent(in) :: i,n
        real(r64), intent(in) :: h
        real(r64) :: f,x
        x = real(n - i, r64) / real(n - 1, r64)
        f = x * h
    end function

    elemental function space_linlog(i,n,h) result(f)
        integer, intent(in) :: i,n
        real(r64), intent(in) :: h
        real(r64) :: f,x
        x = real(i - 1, r64) / real(n - 1, r64)
        f = (1 + h) ** x - 1
    end function

    elemental function space_linlog_r(i,n,h) result(f)
        integer, intent(in) :: i,n
        real(r64), intent(in) :: h
        real(r64) :: f,x
        x = real(n - i, r64) / real(n - 1, r64)
        f = (1 + h) ** x - 1
    end function

    elemental function space_asinh(i,n,h) result(f)
        integer, intent(in) :: i,n
        real(r64), intent(in) :: h
        real(r64) :: f,x
        x = real(i - 1, r64) / real(n - 1, r64)
        f = sinh( x * asinh(h) )
    end function

    elemental function space_asinh_r(i,n,h) result(f)
        integer, intent(in) :: i,n
        real(r64), intent(in) :: h
        real(r64) :: f,x
        x = real(n - i, r64) / real(n - 1, r64)
        f = sinh( x * asinh(h) )
    end function


end module
