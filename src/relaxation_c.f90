module relaxation_c

  use iso_c_binding
  use iso_fortran_env, only: r64 => real64
  use relaxation

  implicit none

contains

  subroutine mrx_number_c(cor, mag, cnd, nr) bind(C, name = 'mrx_model')
    integer(c_int), intent(in), value :: cor, mag, cnd
    integer(c_int), intent(out) :: nr

    nr = mrx_number(cor .ne. 0, mag .ne. 0, cnd .ne. 0)
  end subroutine

  subroutine mrx_get_ny_c(nr, ny) bind(C, name = 'mrx_get_ny')
    integer(c_int), intent(in) :: nr
    integer(c_int), intent(out) :: ny
    integer :: nbl, nbr

    call mrx_sel_dims(nr, ny, nbl, nbr)

  end subroutine

  subroutine mrx_init_c(mb, md, r, alph, zet) bind(C, name = 'mrx_init')
    real(c_double), intent(in), value :: mb, md, r, alph, zet
    mbh = mb
    mdot = md
    radius = r
    alpha = alph
    zeta = zet
  end subroutine


  subroutine mrx_matrix_c(nr,X,nx,Y,ny,A,na,M) &
        & bind(C, name='mrx_matrix')

    integer(c_int), intent(in), value :: nx,ny,na,nr
    real(c_double), intent(in), dimension(nx) :: X
    real(c_double), intent(in), dimension(nx*ny) :: Y
    real(c_double), intent(out), dimension(nx*na) :: A
    real(c_double), intent(out), dimension(nx*na,nx*ny) :: M

    call mrx_matrix(nr,x,Y,M,A)

  end subroutine

  subroutine mrx_advance_c(nr,X,nx,Y,dY,ny) &
        & bind(C, name='mrx_advance')

    integer(c_int), intent(in), value :: nx,ny,nr
    real(c_double), intent(in), dimension(nx) :: X
    real(c_double), intent(in), dimension(nx*ny) :: Y
    real(c_double), intent(out), dimension(nx*ny) :: dY
    real(r64) :: M(nx*ny,nx*ny)

    call mrx_advance(nr,x,Y,M,dY)

  end subroutine


end module
