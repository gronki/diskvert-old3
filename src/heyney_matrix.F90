module heyney_matrix

    use iso_fortran_env
    use slf_cgs
    use relax_coefficients
    use globals

    implicit none

    interface
    subroutine coeff_type(z, y, D, A, AY, AD, ny) bind(c)
        import real64
        integer, intent(in), value :: ny
        real(real64), intent(in) :: z
        real(real64), intent(in), dimension(ny) :: D
        real(real64), intent(in), dimension(ny) :: Y
        real(real64), intent(out), dimension(ny,ny) :: AD
        real(real64), intent(out), dimension(ny,ny) :: AY
        real(real64), intent(out), dimension(ny) :: A
    end subroutine
    end interface

    procedure(coeff_type), pointer :: coeff_proc => COEFF_SS73DYF

contains

    subroutine generate_coefficients(x,nx,Y,ny,M,A) bind(C)

        integer, intent(in), value :: nx,ny
        real(real64), dimension(nx), intent(in) :: x
        ! uklad: [ Y1(1) Y2(1) Y3(1) Y1(2) Y2(2) Y3(2) ... ]
        real(real64), dimension(nx*ny), intent(in) :: Y
        ! uklad: [ A1(1) A2(1) A3(1) A1(2) A2(2) A3(2) ... ]
        real(real64), dimension(nx*ny), intent(out) :: A
        real(real64), dimension(nx*ny,nx*ny), intent(out) :: M
        real(real64), dimension(ny) :: Am, Ym, Dm
        real(real64), dimension(ny,ny) :: MY, MD
        real(real64) :: dx, xm
        integer :: i

        through_space: do i = 1, nx-1
            dx = x(i+1) - x(i)
            xm = (x(i+1) + x(i)) / 2

            Ym = (Y(ilo(i):ihi(i)) + Y(ilo(i+1):ihi(i+1))) / 2
            Dm = (Y(ilo(i+1):ihi(i+1)) - Y(ilo(i):ihi(i))) / dx

            call coeff_proc(xm, Ym, Dm, Am, MY, MD, ny)

            A(ilo(i) : ihi(i)) = Am
            ! Macierz( nr_rownania, nr_niewiadomej )
            M(ilo(i) : ihi(i), ilo(i) : ihi(i)) = MY / 2 - MD / dx
            M(ilo(i) : ihi(i), ilo(i+1) : ihi(i+1)) = MY / 2 + MD / dx
        end do through_space

    contains

        elemental integer function ilo(i)
            integer, intent(in) :: i
            ilo = 1 + ny*(i-1)
        end function
        elemental integer function ihi(i)
            integer, intent(in) :: i
            ihi = ny*i
        end function

    end subroutine


end module heyney_matrix
