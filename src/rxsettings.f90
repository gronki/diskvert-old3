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

      ! write every iteration to another file? useful to demonstrate relaxation
      case ("-write-all","-all")
        cfg_write_all_iters = .TRUE.

      ! recalculate the cooling-heating balance after relaxation? works best
      ! with -compton switch or alone (not much sense with -corona switch)
      case ("-post-corona")
        cfg_post_corona = .TRUE.
      case ("-no-post-corona")
        cfg_post_corona = .FALSE.

      ! include relativictic term in Compton source function? may cause
      ! some inconsistencies.
      case ("-relativistic", "-rel", "-relcompt")
        use_precise_balance = .TRUE.
      case ("-no-relativistic", "-no-rel", "-no-relcompt")
        use_precise_balance = .FALSE.

      ! enable PP condition for MRI shutdown? can cause trouble for convergence
      case ("-quench","-quench-mri","-qmri")
        use_quench_mri = .TRUE.
      case ("-no-quench","-no-quench-mri","-no-qmri")
        use_quench_mri = .FALSE.

      ! use P_rad in alpha prescription?
      case ("-prad-alpha", "-alphaprad")
        use_prad_in_alpha = .TRUE.
      case ("-no-prad-alpha","-no-alphaprad")
        use_prad_in_alpha = .FALSE.

      end select
    end do
  end subroutine

end module rxsettings
