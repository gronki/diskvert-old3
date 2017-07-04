module fileunits

    use iso_fortran_env, only: input_unit, output_unit, error_unit
    implicit none

    integer :: uinp = input_unit
    integer :: uout = output_unit
    integer :: uerr = error_unit
    integer :: ulog = error_unit
    integer :: upar = output_unit

    character(len=*), parameter :: fmpari   = '(A16,1X,I11)'
    character(len=*), parameter :: fmparic  = '(A16,1X,I11,1X,"#",1X,A)'

    character(len=*), parameter :: fmparf   = '(A16,1X,F11.6)'
    character(len=*), parameter :: fmparfc  = '(A16,1X,F11.6,1X,"#",1X,A)'

    character(len=*), parameter :: fmpare   = '(A16,1X,Es11.4)'
    character(len=*), parameter :: fmparec  = '(A16,1X,Es11.4,1X,"#",1X,A)'

    character(len=*), parameter :: fmpars   = '(A16,1X,11A)'
    character(len=*), parameter :: fmparsc  = '(A16,1X,11A,1X,"#",1X,A)'

    character(len=*), parameter :: fmparl   = '(A16,1X,L11)'
    character(len=*), parameter :: fmparlc  = '(A16,1X,L11,1X,"#",1X,A)'

end module fileunits
