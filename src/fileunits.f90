module fileunits

    use iso_fortran_env, only: input_unit, output_unit, error_unit
    implicit none

    integer :: uinp = input_unit
    integer :: uout = output_unit
    integer :: uerr = error_unit
    integer :: ulog = error_unit
    integer :: upar = output_unit

end module fileunits
