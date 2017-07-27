module relaxation_c

  use iso_c_binding
  use iso_fortran_env, only: r64 => real64
  use relaxation

  implicit none

contains

  subroutine mrx_number_c(cor, mag, cnd, nr) bind(C)
    integer(c_int), intent(in), value :: mag, cnd
    character(c_char), intent(in), value :: cor
    integer(c_int), intent(out) :: nr

    nr = mrx_number(cor, mag .ne. 0, cnd .ne. 0)
  end subroutine

  subroutine mrx_get_ny_c(nr, ny) bind(C)
    integer(c_int), intent(in), value :: nr
    integer(c_int), intent(out) :: ny
    integer :: na, nc, nbl, nbr

    call mrx_sel_dims(nr, ny, nc, na, nbl, nbr)

  end subroutine

  subroutine mrx_init_c(mb, md, r, alph, zet) bind(C)
    real(c_double), intent(in), value :: mb, md, r, alph, zet
    mbh = mb
    mdot = md
    radius = r
    alpha = alph
    zeta = zet
    call cylinder(mbh, mdot, radius, omega, facc, teff, zscale)
  end subroutine


  subroutine mrx_matrix_c(nr,X,nx,Y,ny,A,na,M) bind(C)

    integer(c_int), intent(in), value :: nx,ny,na,nr
    real(c_double), intent(in), dimension(nx) :: X
    real(c_double), intent(in), dimension(nx*ny) :: Y
    real(c_double), intent(out), dimension(nx*na) :: A
    real(c_double), intent(out), dimension(nx*na,nx*ny) :: M

    call mrx_matrix(nr,x,Y,M,A)

  end subroutine

  subroutine mrx_advance_c(nr,X,nx,Y,dY,ny,errno) bind(C)

    integer(c_int), intent(in), value :: nx,ny,nr
    integer(c_int), intent(out) :: errno
    real(c_double), intent(in), dimension(nx) :: X
    real(c_double), intent(in), dimension(nx*ny) :: Y
    real(c_double), intent(out), dimension(nx*ny) :: dY
    integer, dimension(nx*ny) :: ipiv
    real(r64) :: M(nx*ny,nx*ny)

    call mrx_matrix(nr, x, Y, M, dY)
    call dgesv(size(M,2), 1, M, size(M,1), ipiv, dY, size(dY), errno)

  end subroutine


end module
