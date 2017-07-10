module rxsettings

  implicit none

  logical :: cfg_write_all_iters = .FALSE.

contains

  subroutine rdargvrx
    integer :: i,errno
    character(2**8) :: arg

    do i = 1, command_argument_count()
      call get_command_argument(i, arg)
      select case (arg)
      case ("-write-all","-all")
        cfg_write_all_iters = .TRUE.
      end select
    end do
  end subroutine

end module rxsettings
