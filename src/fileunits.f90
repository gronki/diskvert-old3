module fileunits

    use iso_fortran_env,        &
    only: uin => input_unit,    &
          uout => output_unit,  &
          uerr => error_unit
    implicit none

    character(*), parameter :: fmpari   = '(A16,1X,I12)'
    character(*), parameter :: fmparic  = '(A16,1X,I12,1X,"#",1X,A)'

    character(*), parameter :: fmparf   = '(A16,1X,F12.6)'
    character(*), parameter :: fmparfc  = '(A16,1X,F12.6,1X,"#",1X,A)'

    character(*), parameter :: fmpare   = '(A16,1X,Es12.4E3)'
    character(*), parameter :: fmparec  = '(A16,1X,Es12.4E3,1X,"#",1X,A)'

    character(*), parameter :: fmpars   = '(A16,1X,12A)'
    character(*), parameter :: fmparsc  = '(A16,1X,12A,1X,"#",1X,A)'

    character(*), parameter :: fmparl   = '(A16,1X,L12)'
    character(*), parameter :: fmparlc  = '(A16,1X,L12,1X,"#",1X,A)'

    character(*), parameter :: fmhdr  = '("#",1X,13("-"),1X,A)'

    character(*) , parameter :: fmcol = '(2A16)'


end module fileunits
