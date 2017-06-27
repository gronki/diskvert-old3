module grid

    use precision

    implicit none

contains

    elemental function space_linear(i,n,h) result(f)
        integer, intent(in) :: i,n
        real(fp), intent(in) :: h
        real(fp) :: f,x
        x = real(i - 1, fp) / real(n - 1, fp)
        f = x * h
    end function

    elemental function space_linear_r(i,n,h) result(f)
        integer, intent(in) :: i,n
        real(fp), intent(in) :: h
        real(fp) :: f,x
        x = real(n - i, fp) / real(n - 1, fp)
        f = x * h
    end function

    elemental function space_linlog(i,n,h) result(f)
        integer, intent(in) :: i,n
        real(fp), intent(in) :: h
        real(fp) :: f,x
        x = real(i - 1, fp) / real(n - 1, fp)
        f = (1 + h) ** x - 1
    end function

    elemental function space_linlog_r(i,n,h) result(f)
        integer, intent(in) :: i,n
        real(fp), intent(in) :: h
        real(fp) :: f,x
        x = real(n - i, fp) / real(n - 1, fp)
        f = (1 + h) ** x - 1
    end function

    elemental function space_asinh(i,n,h) result(f)
        integer, intent(in) :: i,n
        real(fp), intent(in) :: h
        real(fp) :: f,x
        x = real(i - 1, fp) / real(n - 1, fp)
        f = sinh( x * asinh(h) )
    end function

    elemental function space_asinh_r(i,n,h) result(f)
        integer, intent(in) :: i,n
        real(fp), intent(in) :: h
        real(fp) :: f,x
        x = real(n - i, fp) / real(n - 1, fp)
        f = sinh( x * asinh(h) )
    end function


end module
