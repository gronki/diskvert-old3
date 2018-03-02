module rxsettings

  use relaxation
  implicit none

  logical :: cfg_write_all_iters = .FALSE.
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

      case ("-post-corona")
        cfg_post_corona = .TRUE.
      case ("-no-post-corona")
        cfg_post_corona = .FALSE.

      case ("-relativistic", "-rel", "-relcompt")
        use_precise_balance = .TRUE.
      case ("-no-relativistic", "-no-rel", "-no-relcompt")
        use_precise_balance = .FALSE.

      case ("-quench","-quench-mri","-qmri")
        use_quench_mri = .TRUE.
      case ("-no-quench","-no-quench-mri","-no-qmri")
        use_quench_mri = .FALSE.

      end select
    end do
  end subroutine

end module rxsettings
