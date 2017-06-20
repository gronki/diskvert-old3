module setup

    use settings
    use globals
    use confort
    use iso_fortran_env

    implicit none

contains

!------------------------------------ INIT ------------------------------------!
!   reads the configuration, command line arguments, calls the conf and cli    !
!          handlers for specific code, and finally evaluates globals           !
!----------------------------------- INPUTS -----------------------------------!
!                 conf_reader: configuration reading procedure                 !
!                 argv_reader: cli arguments reading procedure                 !
!---------------------------------- OUTPUTS -----------------------------------!
!                              errno: error code                               !
!------------------------------------------------------------------------------!

    subroutine init(errno, conf_reader, argv_reader)

        type(config) :: cfg

        integer, intent(inout) :: errno
        interface
            subroutine cfghandler(cfg,errno)
                import config
                type(config), intent(inout) :: cfg
                integer, intent(inout) :: errno
            end subroutine
            subroutine clhandler(errno)
                integer, intent(inout) :: errno
            end subroutine
        end interface
        procedure(cfghandler), optional :: conf_reader
        procedure(clhandler), optional :: argv_reader
        character(len=150) :: progname

        write(progname, '("DISKVERT v.",I0)') VERSION

        write (error_unit,'(A)')    "",    &
        &   "                     " // trim(progname), &
        &   "     Computes vertical structure of an accretion disk.         ", &
        &   "          Copyright (c) 2017, Dominik Gronkiewicz              ", &
        &   "              Distributed under MIT License.             ", &
        &   "  THIS SOFTWARE IS PROVIDED ""AS IS"", WITHOUT ANY WARRANTY.   ", ""

        call read_command_line(errno)
        if (present(argv_reader)) call argv_reader(errno)
        call mincf_read(cfg)
        if ( cfg % failed() )   error stop "error reading config!"
        call read_globals(cfg)
        if (present(conf_reader)) call conf_reader(cfg,errno)
        call mincf_free(cfg)

        if ( len_trim(outfn_cmdline) .gt. 0 ) then
            outfn = outfn_cmdline
        else if ( len_trim(outfn_config) .gt. 0 ) then
            outfn = outfn_config
        end if

        call eval_globals

    end subroutine

end module
