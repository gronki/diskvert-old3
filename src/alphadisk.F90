module alphadisk

    use iso_fortran_env
    use globals
    use slf_cgs

    implicit none

    integer, parameter, private :: nn = 1024
    integer, parameter :: alphadisk_n = nn

    real(real64), dimension(nn), target, private :: z
    real(real64), dimension(4,nn), target, private :: y, dy
    integer, parameter, private :: c_rho = 1, c_temp = 2, c_flux = 3, c_tau = 4
    real(real64), pointer, dimension(:) :: alphadisk_rho
    real(real64), pointer, dimension(:) :: alphadisk_temp
    real(real64), pointer, dimension(:) :: alphadisk_flux
    real(real64), pointer, dimension(:) :: alphadisk_tau
    real(real64), pointer, dimension(:) :: alphadisk_z

contains

    subroutine quick_alpha_disk(alpha)

        real(real64), intent(in) :: alpha
        real(real64) :: h, h_hi, h_lo
        integer :: it,nmax

        alphadisk_rho  => y(c_rho,nn:1:-1)
        alphadisk_temp => y(c_temp,nn:1:-1)
        alphadisk_flux => y(c_flux,nn:1:-1)
        alphadisk_tau => y(c_tau,nn:1:-1)
        alphadisk_z => z(nn:1:-1)

        h_hi = 12.0 * 31
        h_lo = 12.0 / 31

        write (*,'(A10,A10,A12,A7,A5)') "ITER", "H", "FLX", "MAX", "ACTN"

        adjust_height : do it=1,256
            h = sqrt(h_hi*h_lo)

            call grid(IOR(cfg_grid,GRID_REVERSE),h,z,ngrid)

            y(c_flux,1) = flux_acc
            y(c_rho,1) = 1d-16
            y(c_tau,1) = epsilon(1d0)
            y(c_temp,1) = temp_eff * 0.840896415

            call rk4integr(z, y(:,1), f, y, dy, nmax)

            if ( nmax < nn ) then
                h_hi = h
                write (*,'(I10,F10.3,Es12.4,F7.1,A5)') it, h, y(c_flux,nmax) / flux_acc, 1d2 * nmax / nn, '\/'
            else if (y(c_flux,nn) > 1d-12 * flux_acc) then
                h_lo = h
                write (*,'(I10,F10.3,Es12.4,F7.1,A5)') it, h, y(c_flux,nmax) / flux_acc, 1d2 * nmax / nn, '/\'
            else
                exit adjust_height
            end if
        end do adjust_height
    contains
        subroutine f(z,y,dy,abort)

            real(real64), intent(in) :: z, y(:)
            real(real64), intent(inout) :: dy(size(y))
            real(real64) :: dpgas, pgas, prad
            logical, intent(inout) :: abort

            if ( any(y < 0) ) then
                abort = .TRUE.
            end if

            pgas = 2 * cgs_k_over_mh * y(c_rho) * y(c_temp)
            prad = cgs_a / 3 * y(c_temp)**4

            dy(c_flux) = alpha * omega * ( pgas + prad )
            dpgas = - y(c_rho) * (omega**2 * z - cgs_kapes / cgs_c * y(c_flux))
            dy(c_temp) = - 3 * cgs_kapes * y(c_rho) &
            & / ( 16 * cgs_stef * y(c_temp)**3 ) &
            & * y(c_flux)
            dy(c_rho) = dpgas / (2 * cgs_k_over_mh * y(c_temp)) &
            & - y(c_rho) / y(c_temp) * dy(c_temp)
            dy(c_tau) = - cgs_kapes * y(c_rho)

        end subroutine

    end subroutine


end module
