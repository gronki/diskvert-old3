module slf_cgs

    use iso_fortran_env

    use precision

    implicit none

    real(fp), parameter, private :: pi = 4*atan(real(1,fp))

    real(fp), parameter :: &
            cgs_boltz = 1.3806581212d-16,     &
            cgs_c = 2.99792458d10, &
            cgs_h = 6.62607554040d-27, &
            cgs_graw = 6.67259858585d-8, &
            cgs_mel = 9.10938975454d-28,    &
            cgs_hbar= 0.5 * cgs_h / pi, &
            cgs_stef = 2 * pi**5 * cgs_boltz**4 / ( 15 * cgs_h**3 * cgs_c**2 ), &
            cgs_a = 4 * cgs_stef / cgs_c, &
            cgs_alpha = 1 / 137.036d0, &
            cgs_qel = sqrt( cgs_alpha * cgs_hbar * cgs_c ), &
            cgs_thomson = 8 * pi / 3 * cgs_qel**4 / ( cgs_c**4 * cgs_mel**2 ), &
            cgs_elrad = cgs_qel**2 / ( cgs_c**2 * cgs_mel ), &
            cgs_mhydr = 1.6733d-24,    &
            cgs_kapes1 = cgs_thomson / cgs_mhydr, &
            cgs_k = cgs_boltz,   &
            cgs_k_over_mh = cgs_boltz / cgs_mhydr , &
            cgs_mh_over_k = cgs_mhydr / cgs_boltz, &
            cgs_k_over_mec = cgs_boltz / ( cgs_mel * cgs_c ), &
            cgs_mec_over_k = ( cgs_mel * cgs_c ) / cgs_boltz, &
            cgs_k_over_mec2 = cgs_boltz / ( cgs_mel * cgs_c**2 ), &
            cgs_mec2_over_k = ( cgs_mel * cgs_c**2 ) / cgs_boltz, &
            cgs_msun = 1.99d33, &
            cgs_lsun = 3.9d33,     &
            cgs_kapes = 0.34d0

    real(fp), parameter :: &
            un_keV_in_kelvin = 8.617328149741d-8,   &
            un_kelvin_in_keV = 1 / 8.617328149741d-8,   &
            un_angstrom_x_keV = 12.4

    real(fp), parameter ::  &
        sol_mass = cgs_msun, sol_lum = cgs_lsun,                &
        sol_rschw = 2 * cgs_graw * cgs_msun / cgs_c**2,         &
        sol_mdot_crit = 4 * pi * (cgs_graw * cgs_msun / cgs_c)  &
                * (cgs_mhydr / cgs_thomson),                    &
        sol_mdot_edd = 12 * sol_mdot_crit,                      &
        sol_facc_edd = 3 * cgs_graw * cgs_msun * sol_mdot_edd   &
            / ( 8 * pi * sol_rschw**3 ),                        &
        sol_omega = sqrt( cgs_graw * cgs_msun / sol_rschw**3 )

contains
        !
        ! solves the temperature when pgas+prad
        !
        real(fp) function FINDTEMPERATURE(prs, ro) result(temp)
            real(fp), intent(in) :: prs, ro
            real(fp)  :: temp_hi, temp_lo, prs_comp, miu
            integer :: i
            temp_hi = 1e9
            temp_lo = 0
            miu = 0.5
            temp = 0.5 * (temp_hi + temp_lo)
            do i = 1, 54
                prs_comp = ro / (miu*cgs_mhydr) * cgs_boltz * temp      &
                    +  cgs_a / 3. * temp**4
                if (prs_comp > prs) then
                    temp_hi = temp
                elseif (prs_comp < prs) then
                    temp_lo = temp
                else
                    exit
                endif
                temp = 0.5 * (temp_hi + temp_lo)
            enddo

        end function

end module
