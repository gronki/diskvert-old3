MODULE RELAX_COEFFICIENTS

USE IEEE_ARITHMETIC
USE SLF_CGS
USE PRECISION
USE GLOBALS

IMPLICIT NONE
CONTAINS

PURE SUBROUTINE COEFF_SS73DYF(z, Y, D, A, AY, AD, ny, neq) 
INTEGER, INTENT(in) :: neq ! array dimension
INTEGER, INTENT(in) :: ny ! array dimension
real(fp), INTENT(in) :: z
real(fp), INTENT(in), DIMENSION(0:ny - 1) :: D
real(fp), INTENT(in), DIMENSION(0:ny - 1) :: Y
real(fp), INTENT(out), DIMENSION(0:neq - 1, 0:ny - 1) :: AY
real(fp), INTENT(out), DIMENSION(0:neq - 1) :: A
real(fp), INTENT(out), DIMENSION(0:neq - 1, 0:ny - 1) :: AD
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

PURE SUBROUTINE COEFF_SS73DYF_SIZE(ny, neq, nbl, nbr) 
INTEGER, INTENT(out) :: neq
INTEGER, INTENT(out) :: nbl
INTEGER, INTENT(out) :: ny
INTEGER, INTENT(out) :: nbr
ny = 3
neq = 3
nbl = 1
nbr = 2
END SUBROUTINE COEFF_SS73DYF_SIZE

PURE SUBROUTINE COEFF_SS73DYF_BL(z, Y, BL, MBL, ny, nbl) 
INTEGER, INTENT(in) :: nbl ! array dimension
INTEGER, INTENT(in) :: ny ! array dimension
real(fp), INTENT(in), DIMENSION(0:ny - 1) :: Y
real(fp), INTENT(out), DIMENSION(0:nbl - 1, 0:ny - 1) :: MBL
real(fp), INTENT(out), DIMENSION(0:nbl - 1) :: BL
real(fp), INTENT(in) :: z
BL(0) = Y(2)
MBL(0, 0) = 0
MBL(0, 1) = 0
MBL(0, 2) = 1
END SUBROUTINE COEFF_SS73DYF_BL

PURE SUBROUTINE COEFF_SS73DYF_BR(z, Y, BR, MBR, ny, nbr) 
INTEGER, INTENT(in) :: ny ! array dimension
INTEGER, INTENT(in) :: nbr ! array dimension
real(fp), INTENT(in), DIMENSION(0:ny - 1) :: Y
real(fp), INTENT(out), DIMENSION(0:nbr - 1) :: BR
real(fp), INTENT(out), DIMENSION(0:nbr - 1, 0:ny - 1) :: MBR
real(fp), INTENT(in) :: z
BR(0) = -flux_acc + Y(2)
MBR(0, 0) = 0
MBR(0, 1) = 0
MBR(0, 2) = 1
BR(1) = -2*cgs_stef*Y(1)**4 + Y(2)
MBR(1, 0) = 0
MBR(1, 1) = -8*cgs_stef*Y(1)**3
MBR(1, 2) = 1
END SUBROUTINE COEFF_SS73DYF_BR

PURE SUBROUTINE COEFF_MAGNDYF(z, Y, D, A, AY, AD, ny, neq) 
INTEGER, INTENT(in) :: neq ! array dimension
INTEGER, INTENT(in) :: ny ! array dimension
real(fp), INTENT(in) :: z
real(fp), INTENT(in), DIMENSION(0:ny - 1) :: D
real(fp), INTENT(in), DIMENSION(0:ny - 1) :: Y
real(fp), INTENT(out), DIMENSION(0:neq - 1, 0:ny - 1) :: AY
real(fp), INTENT(out), DIMENSION(0:neq - 1) :: A
real(fp), INTENT(out), DIMENSION(0:neq - 1, 0:ny - 1) :: AD
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
A(2) = Omega*z*zeta*D(3) + D(2)
AY(2, 0) = 0
AD(2, 0) = 0
AY(2, 1) = 0
AD(2, 1) = 0
AY(2, 2) = 0
AD(2, 2) = 1
AY(2, 3) = 0
AD(2, 3) = Omega*z*zeta
A(3) = -Omega*alpha*(cgs_boltz*Y(0)*Y(1)/(cgs_mhydr*miu) + Y(3) + (4.0d0 &
      /3.0d0)*cgs_stef*Y(1)**4/cgs_c) + Omega*z*zeta*D(3) + 2*Omega* &
      zeta*Y(3)
AY(3, 0) = -Omega*alpha*cgs_boltz*Y(1)/(cgs_mhydr*miu)
AD(3, 0) = 0
AY(3, 1) = -Omega*alpha*(cgs_boltz*Y(0)/(cgs_mhydr*miu) + (16.0d0/3.0d0) &
      *cgs_stef*Y(1)**3/cgs_c)
AD(3, 1) = 0
AY(3, 2) = 0
AD(3, 2) = 0
AY(3, 3) = -Omega*alpha + 2*Omega*zeta
AD(3, 3) = Omega*z*zeta
END SUBROUTINE COEFF_MAGNDYF

PURE SUBROUTINE COEFF_MAGNDYF_SIZE(ny, neq, nbl, nbr) 
INTEGER, INTENT(out) :: neq
INTEGER, INTENT(out) :: nbl
INTEGER, INTENT(out) :: ny
INTEGER, INTENT(out) :: nbr
ny = 4
neq = 4
nbl = 2
nbr = 2
END SUBROUTINE COEFF_MAGNDYF_SIZE

PURE SUBROUTINE COEFF_MAGNDYF_BL(z, Y, BL, MBL, ny, nbl) 
INTEGER, INTENT(in) :: nbl ! array dimension
INTEGER, INTENT(in) :: ny ! array dimension
real(fp), INTENT(in), DIMENSION(0:ny - 1) :: Y
real(fp), INTENT(out), DIMENSION(0:nbl - 1, 0:ny - 1) :: MBL
real(fp), INTENT(out), DIMENSION(0:nbl - 1) :: BL
real(fp), INTENT(in) :: z
BL(0) = Y(2)
MBL(0, 0) = 0
MBL(0, 1) = 0
MBL(0, 2) = 1
MBL(0, 3) = 0
BL(1) = -Omega*alpha*(cgs_boltz*Y(0)*Y(1)/(cgs_mhydr*miu) + Y(3) + ( &
      4.0d0/3.0d0)*cgs_stef*Y(1)**4/cgs_c) + 2*Omega*zeta*Y(3)
MBL(1, 0) = -Omega*alpha*cgs_boltz*Y(1)/(cgs_mhydr*miu)
MBL(1, 1) = -Omega*alpha*(cgs_boltz*Y(0)/(cgs_mhydr*miu) + (16.0d0/3.0d0 &
      )*cgs_stef*Y(1)**3/cgs_c)
MBL(1, 2) = 0
MBL(1, 3) = -Omega*alpha + 2*Omega*zeta
END SUBROUTINE COEFF_MAGNDYF_BL

PURE SUBROUTINE COEFF_MAGNDYF_BR(z, Y, BR, MBR, ny, nbr) 
INTEGER, INTENT(in) :: ny ! array dimension
INTEGER, INTENT(in) :: nbr ! array dimension
real(fp), INTENT(in) :: z
real(fp), INTENT(in), DIMENSION(0:ny - 1) :: Y
real(fp), INTENT(out), DIMENSION(0:nbr - 1) :: BR
real(fp), INTENT(out), DIMENSION(0:nbr - 1, 0:ny - 1) :: MBR
BR(0) = 2*Omega*z*zeta*Y(3) - flux_acc + Y(2)
MBR(0, 0) = 0
MBR(0, 1) = 0
MBR(0, 2) = 1
MBR(0, 3) = 2*Omega*z*zeta
BR(1) = -2*cgs_stef*Y(1)**4 + Y(2)
MBR(1, 0) = 0
MBR(1, 1) = -8*cgs_stef*Y(1)**3
MBR(1, 2) = 1
MBR(1, 3) = 0
END SUBROUTINE COEFF_MAGNDYF_BR

PURE SUBROUTINE COEFF_SS73COR(z, Y, D, A, AY, AD, ny, neq) 
INTEGER, INTENT(in) :: neq ! array dimension
INTEGER, INTENT(in) :: ny ! array dimension
real(fp), INTENT(in) :: z
real(fp), INTENT(in), DIMENSION(0:ny - 1) :: D
real(fp), INTENT(in), DIMENSION(0:ny - 1) :: Y
real(fp), INTENT(out), DIMENSION(0:neq - 1, 0:ny - 1) :: AY
real(fp), INTENT(out), DIMENSION(0:neq - 1) :: A
real(fp), INTENT(out), DIMENSION(0:neq - 1, 0:ny - 1) :: AD
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

PURE SUBROUTINE COEFF_SS73COR_SIZE(ny, neq, nbl, nbr) 
INTEGER, INTENT(out) :: neq
INTEGER, INTENT(out) :: nbl
INTEGER, INTENT(out) :: ny
INTEGER, INTENT(out) :: nbr
ny = 4
neq = 4
nbl = 2
nbr = 2
END SUBROUTINE COEFF_SS73COR_SIZE

PURE SUBROUTINE COEFF_SS73COR_BL(z, Y, BL, MBL, ny, nbl) 
INTEGER, INTENT(in) :: nbl ! array dimension
INTEGER, INTENT(in) :: ny ! array dimension
real(fp), INTENT(in), DIMENSION(0:ny - 1) :: Y
real(fp), INTENT(out), DIMENSION(0:nbl - 1, 0:ny - 1) :: MBL
real(fp), INTENT(out), DIMENSION(0:nbl - 1) :: BL
real(fp), INTENT(in) :: z
BL(0) = Y(2)
MBL(0, 0) = 0
MBL(0, 1) = 0
MBL(0, 2) = 1
MBL(0, 3) = 0
BL(1) = Y(1) - Y(3)
MBL(1, 0) = 0
MBL(1, 1) = 1
MBL(1, 2) = 0
MBL(1, 3) = -1
END SUBROUTINE COEFF_SS73COR_BL

PURE SUBROUTINE COEFF_SS73COR_BR(z, Y, BR, MBR, ny, nbr) 
INTEGER, INTENT(in) :: ny ! array dimension
INTEGER, INTENT(in) :: nbr ! array dimension
real(fp), INTENT(in), DIMENSION(0:ny - 1) :: Y
real(fp), INTENT(out), DIMENSION(0:nbr - 1) :: BR
real(fp), INTENT(out), DIMENSION(0:nbr - 1, 0:ny - 1) :: MBR
real(fp), INTENT(in) :: z
BR(0) = -flux_acc + Y(2)
MBR(0, 0) = 0
MBR(0, 1) = 0
MBR(0, 2) = 1
MBR(0, 3) = 0
BR(1) = -2*cgs_stef*Y(3)**4 + Y(2)
MBR(1, 0) = 0
MBR(1, 1) = 0
MBR(1, 2) = 1
MBR(1, 3) = -8*cgs_stef*Y(3)**3
END SUBROUTINE COEFF_SS73COR_BR

PURE SUBROUTINE COEFF_MAGNCOR(z, Y, D, A, AY, AD, ny, neq) 
INTEGER, INTENT(in) :: neq ! array dimension
INTEGER, INTENT(in) :: ny ! array dimension
real(fp), INTENT(in) :: z
real(fp), INTENT(in), DIMENSION(0:ny - 1) :: D
real(fp), INTENT(in), DIMENSION(0:ny - 1) :: Y
real(fp), INTENT(out), DIMENSION(0:neq - 1, 0:ny - 1) :: AY
real(fp), INTENT(out), DIMENSION(0:neq - 1) :: A
real(fp), INTENT(out), DIMENSION(0:neq - 1, 0:ny - 1) :: AD
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
A(2) = Omega*z*zeta*D(4) + D(2)
AY(2, 0) = 0
AD(2, 0) = 0
AY(2, 1) = 0
AD(2, 1) = 0
AY(2, 2) = 0
AD(2, 2) = 1
AY(2, 3) = 0
AD(2, 3) = 0
AY(2, 4) = 0
AD(2, 4) = Omega*z*zeta
A(3) = -Omega*alpha*(cgs_boltz*Y(0)*Y(1)/(cgs_mhydr*miu) + Y(4) + (4.0d0 &
      /3.0d0)*cgs_stef*Y(3)**4/cgs_c) + Omega*z*zeta*D(4) + 2*Omega* &
      zeta*Y(4)
AY(3, 0) = -Omega*alpha*cgs_boltz*Y(1)/(cgs_mhydr*miu)
AD(3, 0) = 0
AY(3, 1) = -Omega*alpha*cgs_boltz*Y(0)/(cgs_mhydr*miu)
AD(3, 1) = 0
AY(3, 2) = 0
AD(3, 2) = 0
AY(3, 3) = -16.0d0/3.0d0*Omega*alpha*cgs_stef*Y(3)**3/cgs_c
AD(3, 3) = 0
AY(3, 4) = -Omega*alpha + 2*Omega*zeta
AD(3, 4) = Omega*z*zeta
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

PURE SUBROUTINE COEFF_MAGNCOR_SIZE(ny, neq, nbl, nbr) 
INTEGER, INTENT(out) :: neq
INTEGER, INTENT(out) :: nbl
INTEGER, INTENT(out) :: ny
INTEGER, INTENT(out) :: nbr
ny = 5
neq = 5
nbl = 3
nbr = 2
END SUBROUTINE COEFF_MAGNCOR_SIZE

PURE SUBROUTINE COEFF_MAGNCOR_BL(z, Y, BL, MBL, ny, nbl) 
INTEGER, INTENT(in) :: nbl ! array dimension
INTEGER, INTENT(in) :: ny ! array dimension
real(fp), INTENT(in), DIMENSION(0:ny - 1) :: Y
real(fp), INTENT(out), DIMENSION(0:nbl - 1, 0:ny - 1) :: MBL
real(fp), INTENT(out), DIMENSION(0:nbl - 1) :: BL
real(fp), INTENT(in) :: z
BL(0) = Y(2)
MBL(0, 0) = 0
MBL(0, 1) = 0
MBL(0, 2) = 1
MBL(0, 3) = 0
MBL(0, 4) = 0
BL(1) = -Omega*alpha*(cgs_boltz*Y(0)*Y(1)/(cgs_mhydr*miu) + Y(4) + ( &
      4.0d0/3.0d0)*cgs_stef*Y(3)**4/cgs_c) + 2*Omega*zeta*Y(4)
MBL(1, 0) = -Omega*alpha*cgs_boltz*Y(1)/(cgs_mhydr*miu)
MBL(1, 1) = -Omega*alpha*cgs_boltz*Y(0)/(cgs_mhydr*miu)
MBL(1, 2) = 0
MBL(1, 3) = -16.0d0/3.0d0*Omega*alpha*cgs_stef*Y(3)**3/cgs_c
MBL(1, 4) = -Omega*alpha + 2*Omega*zeta
BL(2) = Y(1) - Y(3)
MBL(2, 0) = 0
MBL(2, 1) = 1
MBL(2, 2) = 0
MBL(2, 3) = -1
MBL(2, 4) = 0
END SUBROUTINE COEFF_MAGNCOR_BL

PURE SUBROUTINE COEFF_MAGNCOR_BR(z, Y, BR, MBR, ny, nbr) 
INTEGER, INTENT(in) :: ny ! array dimension
INTEGER, INTENT(in) :: nbr ! array dimension
real(fp), INTENT(in) :: z
real(fp), INTENT(in), DIMENSION(0:ny - 1) :: Y
real(fp), INTENT(out), DIMENSION(0:nbr - 1) :: BR
real(fp), INTENT(out), DIMENSION(0:nbr - 1, 0:ny - 1) :: MBR
BR(0) = 2*Omega*z*zeta*Y(4) - flux_acc + Y(2)
MBR(0, 0) = 0
MBR(0, 1) = 0
MBR(0, 2) = 1
MBR(0, 3) = 0
MBR(0, 4) = 2*Omega*z*zeta
BR(1) = -2*cgs_stef*Y(3)**4 + Y(2)
MBR(1, 0) = 0
MBR(1, 1) = 0
MBR(1, 2) = 1
MBR(1, 3) = -8*cgs_stef*Y(3)**3
MBR(1, 4) = 0
END SUBROUTINE COEFF_MAGNCOR_BR

PURE SUBROUTINE COEFF_MAGNCORCND(z, Y, D, A, AY, AD, ny, neq) 
INTEGER, INTENT(in) :: neq ! array dimension
INTEGER, INTENT(in) :: ny ! array dimension
real(fp), INTENT(in) :: z
real(fp), INTENT(in), DIMENSION(0:ny - 1) :: D
real(fp), INTENT(in), DIMENSION(0:ny - 1) :: Y
real(fp), INTENT(out), DIMENSION(0:neq - 1, 0:ny - 1) :: AY
real(fp), INTENT(out), DIMENSION(0:neq - 1) :: A
real(fp), INTENT(out), DIMENSION(0:neq - 1, 0:ny - 1) :: AD
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
A(2) = Omega*z*zeta*D(4) + D(2) + D(5)
AY(2, 0) = 0
AD(2, 0) = 0
AY(2, 1) = 0
AD(2, 1) = 0
AY(2, 2) = 0
AD(2, 2) = 1
AY(2, 3) = 0
AD(2, 3) = 0
AY(2, 4) = 0
AD(2, 4) = Omega*z*zeta
AY(2, 5) = 0
AD(2, 5) = 1
A(3) = -Omega*alpha*(cgs_boltz*Y(0)*Y(1)/(cgs_mhydr*miu) + Y(4) + (4.0d0 &
      /3.0d0)*cgs_stef*Y(3)**4/cgs_c) + Omega*z*zeta*D(4) + 2*Omega* &
      zeta*Y(4)
AY(3, 0) = -Omega*alpha*cgs_boltz*Y(1)/(cgs_mhydr*miu)
AD(3, 0) = 0
AY(3, 1) = -Omega*alpha*cgs_boltz*Y(0)/(cgs_mhydr*miu)
AD(3, 1) = 0
AY(3, 2) = 0
AD(3, 2) = 0
AY(3, 3) = -16.0d0/3.0d0*Omega*alpha*cgs_stef*Y(3)**3/cgs_c
AD(3, 3) = 0
AY(3, 4) = -Omega*alpha + 2*Omega*zeta
AD(3, 4) = Omega*z*zeta
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

PURE SUBROUTINE COEFF_MAGNCORCND_SIZE(ny, neq, nbl, nbr) 
INTEGER, INTENT(out) :: neq
INTEGER, INTENT(out) :: nbl
INTEGER, INTENT(out) :: ny
INTEGER, INTENT(out) :: nbr
ny = 6
neq = 6
nbl = 4
nbr = 2
END SUBROUTINE COEFF_MAGNCORCND_SIZE

PURE SUBROUTINE COEFF_MAGNCORCND_BL(z, Y, BL, MBL, ny, nbl) 
INTEGER, INTENT(in) :: nbl ! array dimension
INTEGER, INTENT(in) :: ny ! array dimension
real(fp), INTENT(in), DIMENSION(0:ny - 1) :: Y
real(fp), INTENT(out), DIMENSION(0:nbl - 1, 0:ny - 1) :: MBL
real(fp), INTENT(out), DIMENSION(0:nbl - 1) :: BL
real(fp), INTENT(in) :: z
BL(0) = Y(2)
MBL(0, 0) = 0
MBL(0, 1) = 0
MBL(0, 2) = 1
MBL(0, 3) = 0
MBL(0, 4) = 0
MBL(0, 5) = 0
BL(1) = -Omega*alpha*(cgs_boltz*Y(0)*Y(1)/(cgs_mhydr*miu) + Y(4) + ( &
      4.0d0/3.0d0)*cgs_stef*Y(3)**4/cgs_c) + 2*Omega*zeta*Y(4)
MBL(1, 0) = -Omega*alpha*cgs_boltz*Y(1)/(cgs_mhydr*miu)
MBL(1, 1) = -Omega*alpha*cgs_boltz*Y(0)/(cgs_mhydr*miu)
MBL(1, 2) = 0
MBL(1, 3) = -16.0d0/3.0d0*Omega*alpha*cgs_stef*Y(3)**3/cgs_c
MBL(1, 4) = -Omega*alpha + 2*Omega*zeta
MBL(1, 5) = 0
BL(2) = Y(1) - Y(3)
MBL(2, 0) = 0
MBL(2, 1) = 1
MBL(2, 2) = 0
MBL(2, 3) = -1
MBL(2, 4) = 0
MBL(2, 5) = 0
BL(3) = Y(5)
MBL(3, 0) = 0
MBL(3, 1) = 0
MBL(3, 2) = 0
MBL(3, 3) = 0
MBL(3, 4) = 0
MBL(3, 5) = 1
END SUBROUTINE COEFF_MAGNCORCND_BL

PURE SUBROUTINE COEFF_MAGNCORCND_BR(z, Y, BR, MBR, ny, nbr) 
INTEGER, INTENT(in) :: ny ! array dimension
INTEGER, INTENT(in) :: nbr ! array dimension
real(fp), INTENT(in) :: z
real(fp), INTENT(in), DIMENSION(0:ny - 1) :: Y
real(fp), INTENT(out), DIMENSION(0:nbr - 1) :: BR
real(fp), INTENT(out), DIMENSION(0:nbr - 1, 0:ny - 1) :: MBR
BR(0) = 2*Omega*z*zeta*Y(4) - flux_acc + Y(2) + Y(5)
MBR(0, 0) = 0
MBR(0, 1) = 0
MBR(0, 2) = 1
MBR(0, 3) = 0
MBR(0, 4) = 2*Omega*z*zeta
MBR(0, 5) = 1
BR(1) = -2*cgs_stef*Y(3)**4 + Y(2)
MBR(1, 0) = 0
MBR(1, 1) = 0
MBR(1, 2) = 1
MBR(1, 3) = -8*cgs_stef*Y(3)**3
MBR(1, 4) = 0
MBR(1, 5) = 0
END SUBROUTINE COEFF_MAGNCORCND_BR


END MODULE
