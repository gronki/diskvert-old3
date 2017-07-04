module alphasimp

    use iso_fortran_env, only: r64 => real64

    use globals
    use slf_cgs
    use grid
    use slf_rk4integr
    use fileunits

    implicit none

    integer, parameter :: nn = 1024

    real(r64), dimension(nn) :: z
    real(r64), dimension(4,nn) :: y, dy
    integer, parameter :: c_rho = 1, c_temp = 2, &
        c_flux = 3, c_tau = 4

contains

    subroutine run_alpha_simple(mbh,mdot,radius,alpha)

        real(r64), intent(in) :: mbh, mdot, radius, alpha
        real(r64) :: omega, Facc, Teff, zscale
        real(r64) :: h, h_hi, h_lo
        integer :: it,nmax,i

        call cylinder(mbh, mdot, radius, omega, Facc, Teff, zscale)

        h_hi = 12.0 * 31
        h_lo = 12.0 / 31

        write (ulog,'(A10,A10,A12,A7,A5)') "ITER", "H", "FLX", "MAX", "ACTN"

        adjust_height : do it=1,256
            h = sqrt(h_hi*h_lo)

            forall (i = 1:nn) z(i) = space_linlog_r(i,nn,h) * zscale

            y(c_flux,1) = Facc
            y(c_rho,1) = 1d-16
            y(c_tau,1) = epsilon(1d0)
            y(c_temp,1) = Teff * 0.5d0 ** (1d0 / 4d0)

            call rk4integr(z, y(:,1), f, y, dy, nmax)

            if ( nmax < nn ) then
                h_hi = h
                write (ulog,'(I10,F10.3,Es12.4,F7.1,A5)') it, h, y(c_flux,nmax) / Facc, 1d2 * nmax / nn, '\/'
            else if (y(c_flux,nn) > 1d-12 * Facc) then
                h_lo = h
                write (ulog,'(I10,F10.3,Es12.4,F7.1,A5)') it, h, y(c_flux,nmax) / Facc, 1d2 * nmax / nn, '/\'
            else
                exit adjust_height
            end if
        end do adjust_height
    contains
        subroutine f(z,y,dy,abort)

            real(r64), intent(in) :: z, y(:)
            real(r64), intent(inout) :: dy(size(y))
            real(r64) :: dpgas, pgas, prad
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
