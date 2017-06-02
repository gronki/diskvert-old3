module heyney_matrix

    use iso_fortran_env
    use slf_cgs
    use relax_coefficients
    use globals

    implicit none

contains

    subroutine generate_coefficients(x,nx,Y,ny,M,A)
        integer, intent(in) :: nx,ny
        real(real64), dimension(nx), intent(in) :: x
        real(real64), dimension(nx*ny), intent(in) :: Y
        real(real64), dimension(nx*ny), intent(out) :: A
        real(real64), dimension(nx*ny,nx*ny), intent(out) :: M
        real(real64), dimension(ny) :: Am, MY, MD
        real(real64) :: dx, xm
        integer :: i

        do concurrent (i = 1:nx-1)
            dx = x(i+1) - x(i)
            xm = (x(i+1) + x(i)) / 2
            
        end do
    end subroutine


end module heyney_matrix
