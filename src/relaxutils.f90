module relaxutils

  use globals
  use iso_fortran_env, only: dp => real64
  implicit none

contains

  subroutine mktau(z, rho, T, tau)
    real(dp), intent(in), dimension(:) :: z, rho,T
    real(dp), intent(out), dimension(:) :: tau
    real(dp) :: dz, rhom, Tm
    integer :: i
    tau(size(tau)) = 1e-12
    do i = size(tau)-1, 1, -1
      rhom = (rho(i+1) + rho(i)) / 2
      Tm = (T(i+1) + T(i)) / 2
      dz = z(i) - z(i+1)
      tau(i) = tau(i+1) - (fksct(rhom,Tm)+fkabs(rhom,Tm)) * rhom * dz
    end do
  end subroutine

  subroutine mktaues(z, rho, T, tau)
    real(dp), intent(in), dimension(:) :: z, rho,T
    real(dp), intent(out), dimension(:) :: tau
    real(dp) :: dz, rhom, Tm
    integer :: i
    tau(size(tau)) = 1e-12
    do i = size(tau)-1, 1, -1
      rhom = (rho(i+1) + rho(i)) / 2
      Tm = (T(i+1) + T(i)) / 2
      dz = z(i) - z(i+1)
      tau(i) = tau(i+1) - kapes0(abuX,abuZ) * rhom * dz
    end do
  end subroutine

end module relaxutils
