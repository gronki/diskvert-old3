module settings

    use iso_fortran_env, only: r64 => real64

    use results
    use confort
    use globals
    use grid

    implicit none

    character(2**8) :: outfn = "disk"
    integer :: ngrid = 2**10

    integer, parameter ::   GRID_LINEAR = 1, &
                        &   GRID_LOG    = 2, &
                        &   GRID_ASINH  = 3
    integer :: tgrid = GRID_LOG
    real(r64) :: htop = 60

contains

!----------------------------- READ_COMMAND_LINE ------------------------------!
!       reads command line arguments and stores them in global variables       !
!---------------------------------- OUTPUTS -----------------------------------!
!                              errno: error code                               !
!------------------------------------------------------------------------------!

    subroutine RDCMDLN(errno)

        integer, intent(inout) :: errno
        integer :: i
        character(128) :: arg, nextarg

        iterate_arguments : do i=1,command_argument_count()

            call get_command_argument(i,arg)

            select case (arg)
            case("-opacity-ff", "-free-free")
                kramers_opacity_ff = .TRUE.
            case("-no-opacity-ff", "-no-free-free")
                kramers_opacity_ff = .FALSE.

            case("-opacity-bf", "-bound-free")
                kramers_opacity_bf = .TRUE.
            case("-no-opacity-bf", "-no-bound-free")
                kramers_opacity_bf = .FALSE.

            case ("-linear","-grid-linear")
                tgrid = GRID_LINEAR
            case ("-log","-grid-log")
                tgrid = GRID_LOG
            case ("-asinh","-grid-asinh")
                tgrid = GRID_ASINH

            case ("-top","-htop","-z-top")
                call get_command_argument(i+1,nextarg)
                read (nextarg,*,iostat = errno) htop
                if ( errno .ne. 0 ) then
                    error stop "top must be followed by an argument " &
                            & // "(the upper boundary of computation interval)"
                end if

            case ("-N","-n","-ngrid")
                call get_command_argument(i+1,nextarg)
                read (nextarg,*,iostat = errno) ngrid
                if ( errno .ne. 0 ) then
                    error stop "ngrid must be followed by an argument " &
                            & // "(number of bins)"
                end if
            case ("-output","-o")
                call get_command_argument(i+1,outfn)
                if ( len_trim(outfn) .eq. 0 ) then
                    error stop "output must be followed by an argument " &
                            & // "(output filename, without extension)"
                end if
            end select

        end do iterate_arguments

    end subroutine


end module
