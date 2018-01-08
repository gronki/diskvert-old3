module rxsettings

  implicit none

  logical :: cfg_write_all_iters = .FALSE.
  character, parameter :: EQUATION_SIMPBALANCE = 'D'

contains

  subroutine rdargvrx
    integer :: i,errno
    character(2**8) :: arg

    do i = 1, command_argument_count()
      call get_command_argument(i, arg)
      select case (arg)
      case ("-write-all","-all")
        cfg_write_all_iters = .TRUE.
      ! case ("-balance-simple","-corona2")
      !   cfg_temperature_method = EQUATION_SIMPBALANCE
      end select
    end do
  end subroutine

end module rxsettings
