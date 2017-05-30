module slf_kramers

    use iso_fortran_env
    use slf_cgs

    implicit none

contains

    elemental function kram_kappa_es(abu_X,abu_Z) result(kes)
        real(real64), intent(in) :: abu_X, abu_Z
        real(real64) :: kes
        kes = 0.2_real64 * (1 + abu_X)
    end function

    elemental function kram_kappa0_bf(abu_X,abu_Z) result(kbf)
        real(real64), intent(in) :: abu_X, abu_Z
        real(real64) :: kbf
        kbf = 4.34e25_real64 * abu_Z * (1 + abu_X)
    end function

    elemental function kram_kappa0_ff(abu_X,abu_Z) result(kff)
        real(real64), intent(in) :: abu_X, abu_Z
        real(real64) :: kff
        kff = 3.68e22_real64 * (1 - abu_Z) * (1 + abu_X)
    end function

    elemental subroutine kram_kappa0(abu_X, abu_Z, kes, kabs)
        real(real64), intent(in) :: abu_X, abu_Z
        real(real64), intent(out) :: kes,kabs

        kes = kram_kappa_es(abu_X,abu_Z)
        kabs = kram_kappa0_ff(abu_X,abu_Z) + kram_kappa0_bf(abu_X,abu_Z)
    end subroutine


end module
