module fileunits

    use iso_fortran_env, only: input_unit, output_unit, error_unit
    implicit none

    integer :: uinp = input_unit
    integer :: uout = output_unit
    integer :: uerr = error_unit
    integer :: ulog = error_unit
    integer :: upar = output_unit

    character(len=*), parameter :: fmpari   = '(A18,2X,I11)'
    character(len=*), parameter :: fmparic  = '(A18,2X,I11,2X,"#",1X,A)'
    character(len=*), parameter :: fmparf   = '(A18,2X,F11.4)'
    character(len=*), parameter :: fmparfc  = '(A18,2X,F11.4,2X,"#",1X,A)'
    character(len=*), parameter :: fmpare   = '(A18,2X,Es11.4)'
    character(len=*), parameter :: fmparec  = '(A18,2X,Es11.4,2X,"#",1X,A)'
    character(len=*), parameter :: fmpars   = '(A18,2X,A)'
    character(len=*), parameter :: fmparsc  = '(A18,2X,A,2X,"#",1X,A)'

end module fileunits
