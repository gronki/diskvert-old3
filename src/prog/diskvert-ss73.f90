program SS73

    use confort

    use precision
    use globals
    use results
    use setup
    use settings
    use model_ss73

    implicit none

    integer :: errno

    character(len=12), dimension(ny) :: y_labels
    character(len=12), dimension(na) :: a_labels
    integer :: nmax
    real(fp), allocatable, dimension(:) :: z
    real(fp), allocatable, dimension(:,:) :: y,dy,a


    y_labels(c_Pgas)    = 'Pgas'
    y_labels(c_Prad)    = 'Prad'
    y_labels(c_Frad)    = 'Frad'
    y_labels(c_tau)     = 'tau_es'
    a_labels(c_rho)     = 'rho'
    a_labels(c_Tgas)    = 'Tgas'
    a_labels(c_Trad)    = 'Trad'
    a_labels(c_zscal)   = 'h'
    a_labels(c_heat)    = 'Gamma'
    a_labels(c_heat_max)= 'Gammamax'
    a_labels(c_compW)   = 'compW'
    a_labels(c_compY)   = 'compY'

    errno = 0
    call init(errno, conf_reader = disk_read_local)

    open(11, file=(trim(outfn) // '.txt'), action='write')
    call write_globals(11)
    call write_settings(11)

    allocate(z(ngrid))
    allocate(y(ny,ngrid), dy(ny,ngrid), a(na,ngrid))

    call run_ss73(z,ngrid,y,dy,a,nmax)

    write (11, fmt_meta_ec) "rho_0", a(c_rho,ngrid), "Central density"
    write (11, fmt_meta_ec) "temp_0", a(c_tgas,ngrid), "Central gas temperature"
    write (11, fmt_meta_fc) "alpha", alpha, "Alpha parameter"
    close(11)

    call write_results_3(z,ngrid,y,dy,y_labels,ny,a,a_labels,na,outfn)

    deallocate( z,y,dy,a )

contains

    subroutine disk_read_local(cfg,errno)
        integer, intent(inout) :: errno
        type(config), intent(inout) :: cfg
        character(len=2048) :: buf

        call mincf_get(cfg, "alpha", buf)
        if ( cfg % not_found() )  then
            error stop "Magnetic alpha-parameter (key: alpha) is REQUIRED!"
        end if
        read (buf,*) alpha

    end subroutine

end program
