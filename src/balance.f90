module energy_balance

    use iso_fortran_env, only: r64 => real64
    use globals
    use slf_threshold

    implicit none

    real(r64) ::   fcool_heat, &
                & fcool_pgas, &
                & fcool_prad, &
                & fcool_trad

contains


    pure subroutine fcool(logTemp,y,dy)

        real(r64), intent(in) :: logTemp
        real(r64), intent(out) :: y
        real(r64), intent(out), optional :: dy
        real(r64) :: xxA2kesPradk_over_mec, xx4sgT4_m_3cPrad, xx4sgT4, xxA
        real(r64) :: rho,kabs,cool
        real(r64) :: Temp

        Temp = exp(logTemp)

        rho = fcool_pgas * miu * cgs_mhydr / ( cgs_boltz * Temp )
        kabs = fkabs(rho,temp)

        xx4sgT4 = 4 * cgs_stef * Temp**4
        xx4sgT4_m_3cPrad = xx4sgT4 - 3 * cgs_c * fcool_Prad
        xxA2kesPradk_over_mec = 12 * kappa_es * fcool_prad * cgs_k_over_mec

        xxA = kabs * xx4sgT4_m_3cPrad + xxA2kesPradk_over_mec * ( Temp - fcool_Trad )

        cool = xxA * rho

        y = fcool_heat - cool

        if ( present(dy) ) then
            xxA = 4 * (xx4sgT4 - 9 * xx4sgT4_m_3cPrad / 8) * kabs &
                & + xxA2kesPradk_over_mec * Temp
            dy =  - ( rho * xxA - cool )
        end if

    end subroutine

    pure subroutine calc_compswitch(pgas, tgas, compsw)
        real(r64), intent(in) :: tgas,pgas
        real(r64), intent(out) ::compsw
        real(r64) :: rho, epsi, taues, compy, kabs

        rho =  pgas * miu * cgs_mhydr / ( cgs_boltz * tgas )

        kabs = fkabs(rho, tgas)
        epsi = kabs / (kabs + kappa_es)
        taues = (1 - epsi) / sqrt(epsi)
        compy = 4 * cgs_boltz * tgas &
            & / ( cgs_mel * cgs_c**2 ) &
            & * max( taues, taues**2 )
        compsw = thrtanh( log(compy) / log(2d0) )

    end subroutine

    pure subroutine calc_compW(heat, pgas, prad, heat_max, compw)
        real(r64), intent(in) :: heat, pgas, prad
        real(r64), intent(out) :: heat_max, compw

        heat_max = 12 * cgs_kapes * prad  &
            & * ( pgas * miu * cgs_mhydr ) &
            & / ( cgs_mel * cgs_c )

        compw = heat / ( heat_max - heat )
    end subroutine

end module
