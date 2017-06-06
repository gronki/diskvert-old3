module precision

#   ifndef PRECISION
#       define PRECISION 14
#   endif

    use iso_fortran_env

implicit none

    integer, parameter :: fp = selected_real_kind(PRECISION)

end module
