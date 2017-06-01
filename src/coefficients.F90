MODULE RELAX_COEFFICIENTS

USE ISO_FORTRAN_ENV
USE ISO_C_BINDING
USE IEEE_ARITHMETIC
USE SLF_CGS
USE GLOBALS

IMPLICIT NONE
CONTAINS

SUBROUTINE COEFF_SS73DYF(z, Y, D, A, AY, AD, alpha) BIND(C)
REAL(real64), INTENT(in) :: alpha
REAL(real64), INTENT(in) :: z
REAL(real64), INTENT(in), DIMENSION(0:2) :: Y
REAL(real64), INTENT(in), DIMENSION(0:2) :: D
REAL(real64), INTENT(out), DIMENSION(0:2) :: A
REAL(real64), INTENT(out), DIMENSION(0:2, 0:2) :: AD
REAL(real64), INTENT(out), DIMENSION(0:2, 0:2) :: AY
A(0) = Omega**2*z*Y(0) + cgs_boltz*D(0)*Y(1)/(cgs_mhydr*miu) + cgs_boltz &
      *D(1)*Y(0)/(cgs_mhydr*miu) - (cgs_kapes + kappa_abs(Y(0), Y(1)))* &
      Y(0)*Y(2)/cgs_c
AY(0, 0) = Omega**2*z + cgs_boltz*D(1)/(cgs_mhydr*miu) - (cgs_kapes + &
      kappa_abs(Y(0), Y(1)))*Y(2)/cgs_c - kappa_abs_1(Y(0), Y(1))*Y(0)* &
      Y(2)/cgs_c
AD(0, 0) = cgs_boltz*Y(1)/(cgs_mhydr*miu)
AY(0, 1) = cgs_boltz*D(0)/(cgs_mhydr*miu) - kappa_abs_2(Y(0), Y(1))*Y(0) &
      *Y(2)/cgs_c
AD(0, 1) = cgs_boltz*Y(0)/(cgs_mhydr*miu)
AY(0, 2) = -(cgs_kapes + kappa_abs(Y(0), Y(1)))*Y(0)/cgs_c
AD(0, 2) = 0
A(1) = 16*cgs_stef*D(1)*Y(1)**3/((3*cgs_kapes + 3*kappa_abs(Y(0), Y(1))) &
      *Y(0)) + Y(2)
AY(1, 0) = -16*cgs_stef*D(1)*Y(1)**3/((3*cgs_kapes + 3*kappa_abs(Y(0), Y &
      (1)))*Y(0)**2) - 48*cgs_stef*kappa_abs_1(Y(0), Y(1))*D(1)*Y(1)**3 &
      /((3*cgs_kapes + 3*kappa_abs(Y(0), Y(1)))**2*Y(0))
AD(1, 0) = 0
AY(1, 1) = 48*cgs_stef*D(1)*Y(1)**2/((3*cgs_kapes + 3*kappa_abs(Y(0), Y( &
      1)))*Y(0)) - 48*cgs_stef*kappa_abs_2(Y(0), Y(1))*D(1)*Y(1)**3/((3 &
      *cgs_kapes + 3*kappa_abs(Y(0), Y(1)))**2*Y(0))
AD(1, 1) = 16*cgs_stef*Y(1)**3/((3*cgs_kapes + 3*kappa_abs(Y(0), Y(1)))* &
      Y(0))
AY(1, 2) = 1
AD(1, 2) = 0
A(2) = -Omega*alpha*(cgs_boltz*Y(0)*Y(1)/(cgs_mhydr*miu) + (4.0d0/3.0d0) &
      *cgs_stef*Y(1)**4/cgs_c) + D(2)
AY(2, 0) = -Omega*alpha*cgs_boltz*Y(1)/(cgs_mhydr*miu)
AD(2, 0) = 0
AY(2, 1) = -Omega*alpha*(cgs_boltz*Y(0)/(cgs_mhydr*miu) + (16.0d0/3.0d0) &
      *cgs_stef*Y(1)**3/cgs_c)
AD(2, 1) = 0
AY(2, 2) = 0
AD(2, 2) = 1
END SUBROUTINE COEFF_SS73DYF

SUBROUTINE COEFF_MAGNDYF(z, Y, D, A, AY, AD, alpha, eta) BIND(C)
REAL(real64), INTENT(in) :: alpha
REAL(real64), INTENT(in) :: z
REAL(real64), INTENT(in) :: eta
REAL(real64), INTENT(in), DIMENSION(0:3) :: D
REAL(real64), INTENT(in), DIMENSION(0:3) :: Y
REAL(real64), INTENT(out), DIMENSION(0:3) :: A
REAL(real64), INTENT(out), DIMENSION(0:3, 0:3) :: AD
REAL(real64), INTENT(out), DIMENSION(0:3, 0:3) :: AY
A(0) = Omega**2*z*Y(0) + cgs_boltz*D(0)*Y(1)/(cgs_mhydr*miu) + cgs_boltz &
      *D(1)*Y(0)/(cgs_mhydr*miu) + D(3) - (cgs_kapes + kappa_abs(Y(0), &
      Y(1)))*Y(0)*Y(2)/cgs_c
AY(0, 0) = Omega**2*z + cgs_boltz*D(1)/(cgs_mhydr*miu) - (cgs_kapes + &
      kappa_abs(Y(0), Y(1)))*Y(2)/cgs_c - kappa_abs_1(Y(0), Y(1))*Y(0)* &
      Y(2)/cgs_c
AD(0, 0) = cgs_boltz*Y(1)/(cgs_mhydr*miu)
AY(0, 1) = cgs_boltz*D(0)/(cgs_mhydr*miu) - kappa_abs_2(Y(0), Y(1))*Y(0) &
      *Y(2)/cgs_c
AD(0, 1) = cgs_boltz*Y(0)/(cgs_mhydr*miu)
AY(0, 2) = -(cgs_kapes + kappa_abs(Y(0), Y(1)))*Y(0)/cgs_c
AD(0, 2) = 0
AY(0, 3) = 0
AD(0, 3) = 1
A(1) = 16*cgs_stef*D(1)*Y(1)**3/((3*cgs_kapes + 3*kappa_abs(Y(0), Y(1))) &
      *Y(0)) + Y(2)
AY(1, 0) = -16*cgs_stef*D(1)*Y(1)**3/((3*cgs_kapes + 3*kappa_abs(Y(0), Y &
      (1)))*Y(0)**2) - 48*cgs_stef*kappa_abs_1(Y(0), Y(1))*D(1)*Y(1)**3 &
      /((3*cgs_kapes + 3*kappa_abs(Y(0), Y(1)))**2*Y(0))
AD(1, 0) = 0
AY(1, 1) = 48*cgs_stef*D(1)*Y(1)**2/((3*cgs_kapes + 3*kappa_abs(Y(0), Y( &
      1)))*Y(0)) - 48*cgs_stef*kappa_abs_2(Y(0), Y(1))*D(1)*Y(1)**3/((3 &
      *cgs_kapes + 3*kappa_abs(Y(0), Y(1)))**2*Y(0))
AD(1, 1) = 16*cgs_stef*Y(1)**3/((3*cgs_kapes + 3*kappa_abs(Y(0), Y(1)))* &
      Y(0))
AY(1, 2) = 1
AD(1, 2) = 0
AY(1, 3) = 0
AD(1, 3) = 0
A(2) = Omega*eta*z*D(3) + D(2)
AY(2, 0) = 0
AD(2, 0) = 0
AY(2, 1) = 0
AD(2, 1) = 0
AY(2, 2) = 0
AD(2, 2) = 1
AY(2, 3) = 0
AD(2, 3) = Omega*eta*z
A(3) = -Omega*alpha*(cgs_boltz*Y(0)*Y(1)/(cgs_mhydr*miu) + Y(3) + (4.0d0 &
      /3.0d0)*cgs_stef*Y(1)**4/cgs_c) + Omega*eta*z*D(3) + 2*Omega*eta* &
      Y(3)
AY(3, 0) = -Omega*alpha*cgs_boltz*Y(1)/(cgs_mhydr*miu)
AD(3, 0) = 0
AY(3, 1) = -Omega*alpha*(cgs_boltz*Y(0)/(cgs_mhydr*miu) + (16.0d0/3.0d0) &
      *cgs_stef*Y(1)**3/cgs_c)
AD(3, 1) = 0
AY(3, 2) = 0
AD(3, 2) = 0
AY(3, 3) = -Omega*alpha + 2*Omega*eta
AD(3, 3) = Omega*eta*z
END SUBROUTINE COEFF_MAGNDYF

SUBROUTINE COEFF_SS73COR(z, Y, D, A, AY, AD, alpha) BIND(C)
REAL(real64), INTENT(in) :: alpha
REAL(real64), INTENT(in) :: z
REAL(real64), INTENT(in), DIMENSION(0:3) :: D
REAL(real64), INTENT(in), DIMENSION(0:3) :: Y
REAL(real64), INTENT(out), DIMENSION(0:3) :: A
REAL(real64), INTENT(out), DIMENSION(0:3, 0:3) :: AD
REAL(real64), INTENT(out), DIMENSION(0:3, 0:3) :: AY
A(0) = Omega**2*z*Y(0) + cgs_boltz*D(0)*Y(1)/(cgs_mhydr*miu) + cgs_boltz &
      *D(1)*Y(0)/(cgs_mhydr*miu) - (cgs_kapes + kappa_abs(Y(0), Y(1)))* &
      Y(0)*Y(2)/cgs_c
AY(0, 0) = Omega**2*z + cgs_boltz*D(1)/(cgs_mhydr*miu) - (cgs_kapes + &
      kappa_abs(Y(0), Y(1)))*Y(2)/cgs_c - kappa_abs_1(Y(0), Y(1))*Y(0)* &
      Y(2)/cgs_c
AD(0, 0) = cgs_boltz*Y(1)/(cgs_mhydr*miu)
AY(0, 1) = cgs_boltz*D(0)/(cgs_mhydr*miu) - kappa_abs_2(Y(0), Y(1))*Y(0) &
      *Y(2)/cgs_c
AD(0, 1) = cgs_boltz*Y(0)/(cgs_mhydr*miu)
AY(0, 2) = -(cgs_kapes + kappa_abs(Y(0), Y(1)))*Y(0)/cgs_c
AD(0, 2) = 0
AY(0, 3) = 0
AD(0, 3) = 0
A(1) = 16*cgs_stef*D(3)*Y(3)**3/((3*cgs_kapes + 3*kappa_abs(Y(0), Y(1))) &
      *Y(0)) + Y(2)
AY(1, 0) = -16*cgs_stef*D(3)*Y(3)**3/((3*cgs_kapes + 3*kappa_abs(Y(0), Y &
      (1)))*Y(0)**2) - 48*cgs_stef*kappa_abs_1(Y(0), Y(1))*D(3)*Y(3)**3 &
      /((3*cgs_kapes + 3*kappa_abs(Y(0), Y(1)))**2*Y(0))
AD(1, 0) = 0
AY(1, 1) = -48*cgs_stef*kappa_abs_2(Y(0), Y(1))*D(3)*Y(3)**3/((3* &
      cgs_kapes + 3*kappa_abs(Y(0), Y(1)))**2*Y(0))
AD(1, 1) = 0
AY(1, 2) = 1
AD(1, 2) = 0
AY(1, 3) = 48*cgs_stef*D(3)*Y(3)**2/((3*cgs_kapes + 3*kappa_abs(Y(0), Y( &
      1)))*Y(0))
AD(1, 3) = 16*cgs_stef*Y(3)**3/((3*cgs_kapes + 3*kappa_abs(Y(0), Y(1)))* &
      Y(0))
A(2) = -Omega*alpha*(cgs_boltz*Y(0)*Y(1)/(cgs_mhydr*miu) + (4.0d0/3.0d0) &
      *cgs_stef*Y(3)**4/cgs_c) + D(2)
AY(2, 0) = -Omega*alpha*cgs_boltz*Y(1)/(cgs_mhydr*miu)
AD(2, 0) = 0
AY(2, 1) = -Omega*alpha*cgs_boltz*Y(0)/(cgs_mhydr*miu)
AD(2, 1) = 0
AY(2, 2) = 0
AD(2, 2) = 1
AY(2, 3) = -16.0d0/3.0d0*Omega*alpha*cgs_stef*Y(3)**3/cgs_c
AD(2, 3) = 0
A(3) = 4*cgs_stef*(4*cgs_boltz*cgs_kapes*Y(3)**4/(cgs_c**2*cgs_mel) + (Y &
      (1) + Y(3))*(Y(1)**2 + Y(3)**2)*kappa_abs(Y(0), Y(1)))*(Y(1) - Y( &
      3))*Y(0) - D(2)
AY(3, 0) = 4*cgs_stef*(4*cgs_boltz*cgs_kapes*Y(3)**4/(cgs_c**2*cgs_mel) &
      + (Y(1) + Y(3))*(Y(1)**2 + Y(3)**2)*kappa_abs(Y(0), Y(1)))*(Y(1) &
      - Y(3)) + 4*cgs_stef*(Y(1) - Y(3))*(Y(1) + Y(3))*(Y(1)**2 + Y(3) &
      **2)*kappa_abs_1(Y(0), Y(1))*Y(0)
AD(3, 0) = 0
AY(3, 1) = 4*cgs_stef*(4*cgs_boltz*cgs_kapes*Y(3)**4/(cgs_c**2*cgs_mel) &
      + (Y(1) + Y(3))*(Y(1)**2 + Y(3)**2)*kappa_abs(Y(0), Y(1)))*Y(0) + &
      4*cgs_stef*(Y(1) - Y(3))*((Y(1) + Y(3))*(Y(1)**2 + Y(3)**2)* &
      kappa_abs_2(Y(0), Y(1)) + 2*(Y(1) + Y(3))*kappa_abs(Y(0), Y(1))*Y &
      (1) + (Y(1)**2 + Y(3)**2)*kappa_abs(Y(0), Y(1)))*Y(0)
AD(3, 1) = 0
AY(3, 2) = 0
AD(3, 2) = -1
AY(3, 3) = -4*cgs_stef*(4*cgs_boltz*cgs_kapes*Y(3)**4/(cgs_c**2*cgs_mel &
      ) + (Y(1) + Y(3))*(Y(1)**2 + Y(3)**2)*kappa_abs(Y(0), Y(1)))*Y(0 &
      ) + 4*cgs_stef*(Y(1) - Y(3))*(16*cgs_boltz*cgs_kapes*Y(3)**3/( &
      cgs_c**2*cgs_mel) + 2*(Y(1) + Y(3))*kappa_abs(Y(0), Y(1))*Y(3) + &
      (Y(1)**2 + Y(3)**2)*kappa_abs(Y(0), Y(1)))*Y(0)
AD(3, 3) = 0
END SUBROUTINE COEFF_SS73COR

SUBROUTINE COEFF_MAGNCOR(z, Y, D, A, AY, AD, alpha, eta) BIND(C)
REAL(real64), INTENT(in) :: alpha
REAL(real64), INTENT(in) :: z
REAL(real64), INTENT(in) :: eta
REAL(real64), INTENT(in), DIMENSION(0:4) :: D
REAL(real64), INTENT(in), DIMENSION(0:4) :: Y
REAL(real64), INTENT(out), DIMENSION(0:4, 0:4) :: AD
REAL(real64), INTENT(out), DIMENSION(0:4) :: A
REAL(real64), INTENT(out), DIMENSION(0:4, 0:4) :: AY
A(0) = Omega**2*z*Y(0) + cgs_boltz*D(0)*Y(1)/(cgs_mhydr*miu) + cgs_boltz &
      *D(1)*Y(0)/(cgs_mhydr*miu) + D(4) - (cgs_kapes + kappa_abs(Y(0), &
      Y(1)))*Y(0)*Y(2)/cgs_c
AY(0, 0) = Omega**2*z + cgs_boltz*D(1)/(cgs_mhydr*miu) - (cgs_kapes + &
      kappa_abs(Y(0), Y(1)))*Y(2)/cgs_c - kappa_abs_1(Y(0), Y(1))*Y(0)* &
      Y(2)/cgs_c
AD(0, 0) = cgs_boltz*Y(1)/(cgs_mhydr*miu)
AY(0, 1) = cgs_boltz*D(0)/(cgs_mhydr*miu) - kappa_abs_2(Y(0), Y(1))*Y(0) &
      *Y(2)/cgs_c
AD(0, 1) = cgs_boltz*Y(0)/(cgs_mhydr*miu)
AY(0, 2) = -(cgs_kapes + kappa_abs(Y(0), Y(1)))*Y(0)/cgs_c
AD(0, 2) = 0
AY(0, 3) = 0
AD(0, 3) = 0
AY(0, 4) = 0
AD(0, 4) = 1
A(1) = 16*cgs_stef*D(3)*Y(3)**3/((3*cgs_kapes + 3*kappa_abs(Y(0), Y(1))) &
      *Y(0)) + Y(2)
AY(1, 0) = -16*cgs_stef*D(3)*Y(3)**3/((3*cgs_kapes + 3*kappa_abs(Y(0), Y &
      (1)))*Y(0)**2) - 48*cgs_stef*kappa_abs_1(Y(0), Y(1))*D(3)*Y(3)**3 &
      /((3*cgs_kapes + 3*kappa_abs(Y(0), Y(1)))**2*Y(0))
AD(1, 0) = 0
AY(1, 1) = -48*cgs_stef*kappa_abs_2(Y(0), Y(1))*D(3)*Y(3)**3/((3* &
      cgs_kapes + 3*kappa_abs(Y(0), Y(1)))**2*Y(0))
AD(1, 1) = 0
AY(1, 2) = 1
AD(1, 2) = 0
AY(1, 3) = 48*cgs_stef*D(3)*Y(3)**2/((3*cgs_kapes + 3*kappa_abs(Y(0), Y( &
      1)))*Y(0))
AD(1, 3) = 16*cgs_stef*Y(3)**3/((3*cgs_kapes + 3*kappa_abs(Y(0), Y(1)))* &
      Y(0))
AY(1, 4) = 0
AD(1, 4) = 0
A(2) = Omega*eta*z*D(4) + D(2)
AY(2, 0) = 0
AD(2, 0) = 0
AY(2, 1) = 0
AD(2, 1) = 0
AY(2, 2) = 0
AD(2, 2) = 1
AY(2, 3) = 0
AD(2, 3) = 0
AY(2, 4) = 0
AD(2, 4) = Omega*eta*z
A(3) = -Omega*alpha*(cgs_boltz*Y(0)*Y(1)/(cgs_mhydr*miu) + Y(4) + (4.0d0 &
      /3.0d0)*cgs_stef*Y(3)**4/cgs_c) + Omega*eta*z*D(4) + 2*Omega*eta* &
      Y(4)
AY(3, 0) = -Omega*alpha*cgs_boltz*Y(1)/(cgs_mhydr*miu)
AD(3, 0) = 0
AY(3, 1) = -Omega*alpha*cgs_boltz*Y(0)/(cgs_mhydr*miu)
AD(3, 1) = 0
AY(3, 2) = 0
AD(3, 2) = 0
AY(3, 3) = -16.0d0/3.0d0*Omega*alpha*cgs_stef*Y(3)**3/cgs_c
AD(3, 3) = 0
AY(3, 4) = -Omega*alpha + 2*Omega*eta
AD(3, 4) = Omega*eta*z
A(4) = 4*cgs_stef*(4*cgs_boltz*cgs_kapes*Y(3)**4/(cgs_c**2*cgs_mel) + (Y &
      (1) + Y(3))*(Y(1)**2 + Y(3)**2)*kappa_abs(Y(0), Y(1)))*(Y(1) - Y( &
      3))*Y(0) - D(2)
AY(4, 0) = 4*cgs_stef*(4*cgs_boltz*cgs_kapes*Y(3)**4/(cgs_c**2*cgs_mel) &
      + (Y(1) + Y(3))*(Y(1)**2 + Y(3)**2)*kappa_abs(Y(0), Y(1)))*(Y(1) &
      - Y(3)) + 4*cgs_stef*(Y(1) - Y(3))*(Y(1) + Y(3))*(Y(1)**2 + Y(3) &
      **2)*kappa_abs_1(Y(0), Y(1))*Y(0)
AD(4, 0) = 0
AY(4, 1) = 4*cgs_stef*(4*cgs_boltz*cgs_kapes*Y(3)**4/(cgs_c**2*cgs_mel) &
      + (Y(1) + Y(3))*(Y(1)**2 + Y(3)**2)*kappa_abs(Y(0), Y(1)))*Y(0) + &
      4*cgs_stef*(Y(1) - Y(3))*((Y(1) + Y(3))*(Y(1)**2 + Y(3)**2)* &
      kappa_abs_2(Y(0), Y(1)) + 2*(Y(1) + Y(3))*kappa_abs(Y(0), Y(1))*Y &
      (1) + (Y(1)**2 + Y(3)**2)*kappa_abs(Y(0), Y(1)))*Y(0)
AD(4, 1) = 0
AY(4, 2) = 0
AD(4, 2) = -1
AY(4, 3) = -4*cgs_stef*(4*cgs_boltz*cgs_kapes*Y(3)**4/(cgs_c**2*cgs_mel &
      ) + (Y(1) + Y(3))*(Y(1)**2 + Y(3)**2)*kappa_abs(Y(0), Y(1)))*Y(0 &
      ) + 4*cgs_stef*(Y(1) - Y(3))*(16*cgs_boltz*cgs_kapes*Y(3)**3/( &
      cgs_c**2*cgs_mel) + 2*(Y(1) + Y(3))*kappa_abs(Y(0), Y(1))*Y(3) + &
      (Y(1)**2 + Y(3)**2)*kappa_abs(Y(0), Y(1)))*Y(0)
AD(4, 3) = 0
AY(4, 4) = 0
AD(4, 4) = 0
END SUBROUTINE COEFF_MAGNCOR

SUBROUTINE COEFF_MAGNCORCND(z, Y, D, A, AY, AD, alpha, eta) BIND(C)
REAL(real64), INTENT(in) :: alpha
REAL(real64), INTENT(in) :: z
REAL(real64), INTENT(in) :: eta
REAL(real64), INTENT(in), DIMENSION(0:5) :: D
REAL(real64), INTENT(in), DIMENSION(0:5) :: Y
REAL(real64), INTENT(out), DIMENSION(0:5, 0:5) :: AD
REAL(real64), INTENT(out), DIMENSION(0:5, 0:5) :: AY
REAL(real64), INTENT(out), DIMENSION(0:5) :: A
A(0) = Omega**2*z*Y(0) + cgs_boltz*D(0)*Y(1)/(cgs_mhydr*miu) + cgs_boltz &
      *D(1)*Y(0)/(cgs_mhydr*miu) + D(4) - (cgs_kapes + kappa_abs(Y(0), &
      Y(1)))*Y(0)*Y(2)/cgs_c
AY(0, 0) = Omega**2*z + cgs_boltz*D(1)/(cgs_mhydr*miu) - (cgs_kapes + &
      kappa_abs(Y(0), Y(1)))*Y(2)/cgs_c - kappa_abs_1(Y(0), Y(1))*Y(0)* &
      Y(2)/cgs_c
AD(0, 0) = cgs_boltz*Y(1)/(cgs_mhydr*miu)
AY(0, 1) = cgs_boltz*D(0)/(cgs_mhydr*miu) - kappa_abs_2(Y(0), Y(1))*Y(0) &
      *Y(2)/cgs_c
AD(0, 1) = cgs_boltz*Y(0)/(cgs_mhydr*miu)
AY(0, 2) = -(cgs_kapes + kappa_abs(Y(0), Y(1)))*Y(0)/cgs_c
AD(0, 2) = 0
AY(0, 3) = 0
AD(0, 3) = 0
AY(0, 4) = 0
AD(0, 4) = 1
AY(0, 5) = 0
AD(0, 5) = 0
A(1) = 16*cgs_stef*D(3)*Y(3)**3/((3*cgs_kapes + 3*kappa_abs(Y(0), Y(1))) &
      *Y(0)) + Y(2)
AY(1, 0) = -16*cgs_stef*D(3)*Y(3)**3/((3*cgs_kapes + 3*kappa_abs(Y(0), Y &
      (1)))*Y(0)**2) - 48*cgs_stef*kappa_abs_1(Y(0), Y(1))*D(3)*Y(3)**3 &
      /((3*cgs_kapes + 3*kappa_abs(Y(0), Y(1)))**2*Y(0))
AD(1, 0) = 0
AY(1, 1) = -48*cgs_stef*kappa_abs_2(Y(0), Y(1))*D(3)*Y(3)**3/((3* &
      cgs_kapes + 3*kappa_abs(Y(0), Y(1)))**2*Y(0))
AD(1, 1) = 0
AY(1, 2) = 1
AD(1, 2) = 0
AY(1, 3) = 48*cgs_stef*D(3)*Y(3)**2/((3*cgs_kapes + 3*kappa_abs(Y(0), Y( &
      1)))*Y(0))
AD(1, 3) = 16*cgs_stef*Y(3)**3/((3*cgs_kapes + 3*kappa_abs(Y(0), Y(1)))* &
      Y(0))
AY(1, 4) = 0
AD(1, 4) = 0
AY(1, 5) = 0
AD(1, 5) = 0
A(2) = Omega*eta*z*D(4) + D(2) + D(5)
AY(2, 0) = 0
AD(2, 0) = 0
AY(2, 1) = 0
AD(2, 1) = 0
AY(2, 2) = 0
AD(2, 2) = 1
AY(2, 3) = 0
AD(2, 3) = 0
AY(2, 4) = 0
AD(2, 4) = Omega*eta*z
AY(2, 5) = 0
AD(2, 5) = 1
A(3) = -Omega*alpha*(cgs_boltz*Y(0)*Y(1)/(cgs_mhydr*miu) + Y(4) + (4.0d0 &
      /3.0d0)*cgs_stef*Y(3)**4/cgs_c) + Omega*eta*z*D(4) + 2*Omega*eta* &
      Y(4)
AY(3, 0) = -Omega*alpha*cgs_boltz*Y(1)/(cgs_mhydr*miu)
AD(3, 0) = 0
AY(3, 1) = -Omega*alpha*cgs_boltz*Y(0)/(cgs_mhydr*miu)
AD(3, 1) = 0
AY(3, 2) = 0
AD(3, 2) = 0
AY(3, 3) = -16.0d0/3.0d0*Omega*alpha*cgs_stef*Y(3)**3/cgs_c
AD(3, 3) = 0
AY(3, 4) = -Omega*alpha + 2*Omega*eta
AD(3, 4) = Omega*eta*z
AY(3, 5) = 0
AD(3, 5) = 0
A(4) = 4*cgs_stef*(4*cgs_boltz*cgs_kapes*Y(3)**4/(cgs_c**2*cgs_mel) + (Y &
      (1) + Y(3))*(Y(1)**2 + Y(3)**2)*kappa_abs(Y(0), Y(1)))*(Y(1) - Y( &
      3))*Y(0) - D(2)
AY(4, 0) = 4*cgs_stef*(4*cgs_boltz*cgs_kapes*Y(3)**4/(cgs_c**2*cgs_mel) &
      + (Y(1) + Y(3))*(Y(1)**2 + Y(3)**2)*kappa_abs(Y(0), Y(1)))*(Y(1) &
      - Y(3)) + 4*cgs_stef*(Y(1) - Y(3))*(Y(1) + Y(3))*(Y(1)**2 + Y(3) &
      **2)*kappa_abs_1(Y(0), Y(1))*Y(0)
AD(4, 0) = 0
AY(4, 1) = 4*cgs_stef*(4*cgs_boltz*cgs_kapes*Y(3)**4/(cgs_c**2*cgs_mel) &
      + (Y(1) + Y(3))*(Y(1)**2 + Y(3)**2)*kappa_abs(Y(0), Y(1)))*Y(0) + &
      4*cgs_stef*(Y(1) - Y(3))*((Y(1) + Y(3))*(Y(1)**2 + Y(3)**2)* &
      kappa_abs_2(Y(0), Y(1)) + 2*(Y(1) + Y(3))*kappa_abs(Y(0), Y(1))*Y &
      (1) + (Y(1)**2 + Y(3)**2)*kappa_abs(Y(0), Y(1)))*Y(0)
AD(4, 1) = 0
AY(4, 2) = 0
AD(4, 2) = -1
AY(4, 3) = -4*cgs_stef*(4*cgs_boltz*cgs_kapes*Y(3)**4/(cgs_c**2*cgs_mel &
      ) + (Y(1) + Y(3))*(Y(1)**2 + Y(3)**2)*kappa_abs(Y(0), Y(1)))*Y(0 &
      ) + 4*cgs_stef*(Y(1) - Y(3))*(16*cgs_boltz*cgs_kapes*Y(3)**3/( &
      cgs_c**2*cgs_mel) + 2*(Y(1) + Y(3))*kappa_abs(Y(0), Y(1))*Y(3) + &
      (Y(1)**2 + Y(3)**2)*kappa_abs(Y(0), Y(1)))*Y(0)
AD(4, 3) = 0
AY(4, 4) = 0
AD(4, 4) = 0
AY(4, 5) = 0
AD(4, 5) = 0
A(5) = kappa_cond(Y(0), Y(1))*D(1) + Y(5)
AY(5, 0) = kappa_cond_1(Y(0), Y(1))*D(1)
AD(5, 0) = 0
AY(5, 1) = kappa_cond_2(Y(0), Y(1))*D(1)
AD(5, 1) = kappa_cond(Y(0), Y(1))
AY(5, 2) = 0
AD(5, 2) = 0
AY(5, 3) = 0
AD(5, 3) = 0
AY(5, 4) = 0
AD(5, 4) = 0
AY(5, 5) = 1
AD(5, 5) = 0
END SUBROUTINE COEFF_MAGNCORCND


END MODULE
