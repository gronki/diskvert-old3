program disk_properties

    use confort

    use globals
    use setup
    use settings

    implicit none

    integer :: errno

    errno = 0
    call init(errno)

    if ( outfn_cmdline /= "" ) then
        open(33, file=trim(outfn_cmdline), action="write")
        call write_globals(33)
        close(33)
    else
        call write_globals(6)
    end if


end program
