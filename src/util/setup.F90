module setup

    use settings
    use globals
    use confort

    implicit none

contains

    subroutine init(errno, conf_reader, argv_reader)

        type (confort_c) :: cfg

        integer, intent(inout) :: errno
        interface
            subroutine cfghandler(cfg,errno)
                import confort_c
                type(confort_c), intent(inout) :: cfg
                integer, intent(inout) :: errno
            end subroutine
            subroutine clhandler(errno)
                integer, intent(inout) :: errno
            end subroutine
        end interface
        procedure(cfghandler), optional :: conf_reader
        procedure(clhandler), optional :: argv_reader

        write (error_unit,'(A)')    "",    &
        &   "                        DISKVERT ", "", &
        &   "     Computes vertical structure of an accretion disk.         ", &
        &   "          Copyright (c) 2017, Dominik Gronkiewicz              ", &
        &   "              Distributed under MIT License.             ", &
        &   "  THIS SOFTWARE IS PROVIDED ""AS IS"", WITHOUT ANY WARRANTY.   ", ""

        call read_command_line(errno)
        if (present(argv_reader)) call argv_reader(errno)
        call mincf_read(cfg, errno)
        if ( errno .ne. MINCF_OK )   error stop "error reading config!"
        call read_globals(cfg,errno)
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
