module slf_integrate

  use iso_fortran_env, only: r32 => real32, r64 => real64
  implicit none

  interface integrate
    module procedure integrate32
    module procedure integrate64
  end interface integrate

contains

  pure real(r32) function integrate32(y,x) result(intg)
    real(r32), dimension(:), intent(in) :: x, y
    integer :: n

    if (size(x) /= size(y)) error stop "size(x) /= size(y)"

    n = size(x)
    intg = sum((y(2:n) + y(1:n-1)) * (x(2:n) - x(1:n-1))) / 2
  end function

  pure real(r64) function integrate64(y,x) result(intg)
    real(r64), dimension(:), intent(in) :: x, y
    integer :: n

    if (size(x) /= size(y)) error stop "size(x) /= size(y)"

    n = size(x)
    intg = sum((y(2:n) + y(1:n-1)) * (x(2:n) - x(1:n-1))) / 2
  end function

end module slf_integrate
