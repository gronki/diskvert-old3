module rxsettings

  implicit none

  logical :: cfg_write_all_iters = .FALSE.
  logical :: cfg_adjust_height_beta = .TRUE.
  character, parameter :: EQUATION_SIMPBALANCE = 'D'
  logical :: cfg_post_corona = .false.

contains

  subroutine rdargvrx
    integer :: i,errno
    character(2**8) :: arg

    do i = 1, command_argument_count()
      call get_command_argument(i, arg)
      select case (arg)
      case ("-write-all","-all")
        cfg_write_all_iters = .TRUE.

      case ("-adjust-height","-adjust")
        cfg_adjust_height_beta = .TRUE.
      case ("-no-adjust-height","-no-adjust")
        cfg_adjust_height_beta = .FALSE.

      case ("-post-corona")
        cfg_post_corona = .TRUE.
      case ("-no-post-corona")
        cfg_post_corona = .FALSE.

      ! case ("-balance-simple","-corona2")
      !   cfg_temperature_method = EQUATION_SIMPBALANCE

      end select
    end do
  end subroutine

end module rxsettings
