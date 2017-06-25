module settings

    use results
    use confort
    use iso_fortran_env, only: input_unit, output_unit, error_unit
    use globals

contains

!----------------------------- COMMAND_LINE_FLAG ------------------------------!
! scans command line arguments for a given flag (such as -foo) or negation of  !
!                     it (such as -no-foo or -disable-foo)                     !
!----------------------------------- INPUTS -----------------------------------!
!              arg: command line arg (this should be looped over)              !
!             flag: name of the flag (eg. "foo")                               !
!---------------------------------- OUTPUTS -----------------------------------!
!     logi: is set to .TRUE. if flag is set or .FALSE. if negation is set      !
!------------------------------------------------------------------------------!

    subroutine command_line_flag(arg,flag,logi)
        character(len=*), intent(in) :: arg,flag
        logical, intent(inout) :: logi
        character(len=128) :: fullflag

        fullflag = "-" // adjustl(flag)

        if ( arg .eq. fullflag  .or. arg .eq. "-enable" // fullflag ) then
            logi = .true.
#           if VERBOSE
                write (error_unit,"('Flag ',A,' is enabled.')") trim(flag)
#           endif
        end if
        if ( arg .eq. "-no" // fullflag  .or. arg .eq. "-disable" // fullflag ) then
            logi = .false.
#           if VERBOSE
                write (error_unit,"('Flag ',A,' is disabled.')") trim(flag)
#           endif
        end if
    end subroutine

!----------------------------- READ_COMMAND_LINE ------------------------------!
!       reads command line arguments and stores them in global variables       !
!---------------------------------- OUTPUTS -----------------------------------!
!                              errno: error code                               !
!------------------------------------------------------------------------------!

    subroutine read_command_line(errno)

        integer, intent(inout) :: errno
        integer :: i
        character(len=128) :: arg, nextarg

        iterate_arguments : do i=1,command_argument_count()

            call get_command_argument(i,arg)

            call command_line_flag(arg,"opacity-ff",cfg_opacity_ff)
            call command_line_flag(arg,"opacity-bf",cfg_opacity_bf)
            call command_line_flag(arg,"dynamo-shutdown",cfg_allow_mri_shutdown)
            call command_line_flag(arg,"single",cfg_single_run)
            call command_line_flag(arg,"compton-term",cfg_compton_term)
            call command_line_flag(arg,"euler-method",cfg_euler_integration)
            call command_line_flag(arg,"conduction",cfg_conduction)

            select case (arg)
            case ("-compton")
                cfg_temperature_method = EQUATION_COMPTON
            case ("-balance")
                cfg_temperature_method = EQUATION_BALANCE
                cfg_balance_multi = .false.
            case ("-balance-multi")
                cfg_temperature_method = EQUATION_BALANCE
                cfg_balance_multi = .true.
            case ("-equilibrium")
                cfg_temperature_method = EQUATION_EQUILIBR
            case ("-linear","-grid-linear")
                cfg_grid = GRID_LINEAR
            case ("-log","-grid-log")
                cfg_grid = GRID_LOG
            case ("-asinh","-grid-asinh")
                cfg_grid = GRID_ASINH
            case ("-max-iteration-error","-precision")
                call get_command_argument(i+1,nextarg)
                if ( len_trim(nextarg) .gt. 0 ) then
                    read (nextarg,*) max_iteration_error
                else
                    error stop "-precision must be followed by an argument" &
                                & // " (relative error to be reached until " &
                                & // " iteration stops)"
                end if
            case ("-top","-htop","-z-top")
                call get_command_argument(i+1,nextarg)
                if ( len_trim(nextarg) .gt. 0 ) then
                    read (nextarg,*) htop
                else
                    error stop "top must be followed by an argument " &
                            & // "(the upper boundary of computation interval)"
                end if
            case ("-N","-n","-ngrid")
                call get_command_argument(i+1,nextarg)
                if ( len_trim(nextarg) .gt. 0 ) then
                    read (nextarg,*) ngrid
                else
                    error stop "ngrid must be followed by an argument " &
                            & // "(number of bins)"
                end if
            case ("-output","-o")
                call get_command_argument(i+1,nextarg)
                if ( len_trim(nextarg) .gt. 0 ) then
                    read (nextarg,*) outfn_cmdline
                else
                    error stop "output must be followed by an argument " &
                            & // "(output filename, without extension)"
                end if
            end select

        end do iterate_arguments

    end subroutine read_command_line

!------------------------------- WRITE_SETTINGS -------------------------------!
!                  writes global settings to given file unit                   !
!----------------------------------- INPUTS -----------------------------------!
!                        lun: file unit (must be open!)                        !
!------------------------------------------------------------------------------!

    subroutine write_settings(lun)
        integer, intent(in) :: lun

        write (lun, fmt_meta_ic) "N", ngrid, "Liczba elementow"

        select case (cfg_temperature_method)
        case (EQUATION_EQUILIBR)
            write (lun, fmt_meta_s) "equation", "equilibrium"
            write (lun, fmt_meta_s) "equation_description", &
            & "Using radiation temperature"
        case (EQUATION_COMPTON)
            write (lun, fmt_meta_s) "equation", "compton"
            write (lun, fmt_meta_s) "equation_description", &
            & "Assuming compton cooling when Y>1"
        case (EQUATION_BALANCE)
            write (lun, fmt_meta_s) "equation", "balance"
            write (lun, fmt_meta_s) "equation_description", &
            & "Solving heating-cooling balance"
        case default
            error stop "incorrect value of cfg_temperature_method"
        end select

        select case (cfg_grid)
        case (GRID_LINEAR)
            write (lun, fmt_meta_s) "grid_type", "linear"
        case (GRID_LOG)
            write (lun, fmt_meta_s) "grid_type", "log"
        case (GRID_ASINH)
            write (lun, fmt_meta_s) "grid_type", "asinh"
        case default
            error stop "incorrect value of cfg_grid"
        end select
    contains
        subroutine write_flag
        end subroutine
    end subroutine

!-------------------------------- READ_GLOBALS --------------------------------!
!        reads neccesary global parameters from the configuration file         !
!----------------------------------- INPUTS -----------------------------------!
!                          cfg: configuration object                           !
!------------------------------------------------------------------------------!

    subroutine read_globals(cfg)

        type(config), intent(inout) :: cfg
        character(len=2048) :: buf

        errno = 0

        call mincf_get(cfg,'output',outfn_config)

        call mincf_get(cfg,'mass',buf)
        if ( cfg % not_found() ) then
            error stop "Black hole mass (key: mass) is REQUIRED!"
        end if
        read (buf,*) m_bh

        call mincf_get(cfg,'radius',buf)
        if ( cfg % not_found() ) then
            error stop "Distance from BH in Schwarzchild radii " &
                    & // "(key: radius) is REQUIRED!"
        end if
        read (buf,*) r_calc

        call mincf_get(cfg,'accretion_rate',buf)
        if ( cfg % not_found() ) then
            error stop "Accretion rate " &
                    & // "(key: accretion_rate) is REQUIRED!"
        end if
        read (buf,*) m_dot

        call MINCF_GET(cfg, 'abun_Z', buf, '0.02')
        read (buf,*) abun_Z

    end subroutine
    
!------------------------------- WRITE_GLOBALS --------------------------------!
!               writes global values to a given (open) file unit               !
!----------------------------------- INPUTS -----------------------------------!
!                                lun: file unit                                !
!------------------------------------------------------------------------------!

    subroutine write_globals(lun)
        integer, intent(in) :: lun

        write (lun, fmt_meta_ec) "mass",(m_bh),"x Msun"
        write (lun, fmt_meta_ec) "accretion_rate",(m_dot),"x medd"
        write (lun, fmt_meta_ec) "radius_rschw", r_calc,"x Rschw"
        write (lun, fmt_meta_ec) "radius", r_calc*rschw, "cm"

        write (lun, fmt_meta_ec) "rschw", rschw, "cm"
        write (lun, fmt_meta_ec) "rgraw", rgraw, "cm"
        write (lun, fmt_meta_ec) "rschw_km", (rschw / 1e5), "km"

        write (lun, fmt_meta_fc) "rotation_velocity", &
                & (omega*rschw*r_calc / cgs_c), "c"

        write (lun, fmt_meta_ec) "omega", omega, "Angular velocity"
        write (lun, fmt_meta_ec) "csound", csound, "Sound speed"

        write (lun, fmt_meta_ec) "temp_eff", temp_eff, "Effective temperature"
        write (lun, fmt_meta_ec) "flux_acc", flux_acc, "erg/cm2/s"
        write (lun, fmt_meta_e) "kram_abs_0", kram_abs
        write (lun, fmt_meta_e) "kram_es", kram_es
        write (lun, fmt_meta_e) "kram_ff_0", kram_ff
        write (lun, fmt_meta_e) "kram_bf_0", kram_bf
        write (lun, fmt_meta_ec) "zscale", zscale, "cm"
        write (lun, fmt_meta_ec) "zscale_km", (zscale / 1e5), "km"
        write (lun, fmt_meta_fc) "d", (zscale / (r_calc*rschw)), "disk ratio"
    end subroutine

end module
