module slf_ramps

  use iso_fortran_env, only: r64 => real64
  implicit none

contains

  !----------------------------------------------------------------------------!
  ! this is a linear ramp

  elemental function ramp1(i,n,y0) result(y)
    integer, intent(in) :: i,n
    real(r64), intent(in), optional :: y0
    real(r64) :: y
    y = merge(real(i) / n, 1.0, i .le. n)
    if (present(y0)) y = y0 + (1 - y0) * y
  end function

  !----------------------------------------------------------------------------!
  ! this is a very steep ramp, for later iterations

  elemental function ramp2(i,n,y0) result(y)
    integer, intent(in) :: i,n
    real(r64), intent(in), optional :: y0
    real(r64) :: t,y
    t = merge(real(i) / n, 1.0, i .le. n)
    y = t * (2 - t)
    if (present(y0)) y = y0 + (1 - y0) * y
  end function

  !----------------------------------------------------------------------------!
  ! this is a very smooth, s-shaped ramp
  ! slow convergence but good if we are away from the solution

  elemental function ramp3(i,n,y0) result(y)
    integer, intent(in) :: i,n
    real(r64), intent(in), optional :: y0
    real(r64) :: t,y
    t = merge(real(i) / n, 1.0, i .le. n)
    y = (3 - 2*t) * t**2
    if (present(y0)) y = y0 + (1 - y0) * y
  end function

  !----------------------------------------------------------------------------!
  ! this is a very smooth ramp where you can control the starting
  ! level as well as slope of the starting point. It may be applied
  ! for very unstable equations (A = 0, B = 0) or for a very quick
  ! convergence (A = 0.25, B = 1)

  elemental function ramp4(i,n,A,B) result(y)
    integer, intent(in) :: i,n
    real(r64), intent(in) :: A,B
    real(r64) :: x,y
    x = merge(real(i) / n, 1.0, i .le. n)
    y = A + (1 - A) * B * x           &
      + 3 * (A - 1) * (B - 2) * x**2  &
      - (A - 1) * (3*B - 8) * x**3    &
      + (A - 1) * (B - 3) * x**4
  end function

  !----------------------------------------------------------------------------!
  ! this ramp starts smoothly, rises rapidly, and then slowly saturates to 1

  elemental function ramp5(i,n) result(y)
    integer, intent(in) :: i,n
    real(r64) :: t,y
    t = merge(real(i) / n, 1.0, i .le. n)
    y = (3*t**2 - 8*t + 6) * t**2
  end function

  !----------------------------------------------------------------------------!

end module slf_ramps
