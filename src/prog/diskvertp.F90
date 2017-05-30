program B15

    use ieee_arithmetic
    use slf_cgs
    use slf_threshold
    use confort
    use setup
    use settings_io
    use globals

    implicit none

    integer :: iter_max = 5
    integer :: i,iter=1, last_i

    integer, parameter :: n_vars = 4, v_ro=1, v_temp=4, v_pmag=2, v_flux = 3
    integer, parameter :: n_params = 1, p_mri=1

    real (dp) :: tmp_a(n_vars), tmp_b(n_vars), tmp_c(n_vars), tmp_d(n_vars)

    character (len=*) , parameter :: fmt33e = "(A30,Es15.4,1X,A8)"
    character (len=*) , parameter :: fmt33s = "(A30,A15,1X,A8)"
    character (len=*) , parameter :: fmt33f = "(A30,F15.4,1X,A8)"
    character (len=*) , parameter :: fmt33i = "(A30,I15,1X,A8)"

    real(real64), parameter :: pi = 4*atan(real(1,8))
    character(len=200) :: outfn_ext, outfn_info_ext


    real (dp), allocatable  :: z(:), dz(:), prs(:),  &
            tau(:),   alfven(:), csound1(:), beta(:),  val(:,:), drv(:,:), alfa_eff(:)



    real (dp) ::  z_top_0, z_bottom,  beta_0,  ro_0,   &
             alpha, dzeta,  mri_active=0.,    &
            log_step_precision = 3.6, gamm = 1., poly_k,  mri_tolerance, &
            z_top_start

    character (len=*) , parameter :: fmt1 = "(50es15.5)"
    character (len=*) , parameter :: fmt1s = "(50a15)"

    character (len=1024) :: buf

    integer :: errno
    type(confort_c) :: cfg

    last_i = ngrid


    errno = 0
    call init(errno, conf_reader = disk_read_local)

    outfn_ext = trim(outfn) // '.txt'
    outfn_info_ext = trim(outfn) // '_info.txt'
    open(unit=33, file=outfn_info_ext, action='write')

    dzeta = 0.5 * alpha * ( 1. + beta_0 )

    write (33,'(A,Es10.3,A,F8.2,A,F8.5,A)')  'Mbh = ', m_bh, ' x Msun; R = ', r_calc, &
        ' x Rschw; m* = ', acc_rate, ' x m*edd'
    if ( m_bh > 1000. ) then
        write (33,fmt33e)  'Mass black hole', m_bh, 'x M_sun'
    else
        write (33,fmt33f)  'Mass black hole', m_bh, 'x M_sun'
    endif
    write (33,fmt33e)  'Accretion rate (x m_edd)', acc_rate
    write (33,fmt33f)  'Calculation radius', r_calc, '(x Rschw)'
    write (33,fmt33e)  'Equatorial density', ro_0, '[g/cm3]'


    mri_tolerance = 0.25

    if ( rschw < 1e8 ) then
        write (33,fmt33f)  'Schwarz. radius', rschw / 1e5, 'km'
    else
        write (33,fmt33e)  'Schwarz. radius', rschw / 1e5, 'km'
    endif

    write (33,fmt33e)  'Angular veloc.', 24*3600.*omega/pi, '(1/day)'
    write (33,fmt33f)  'Linear veloc.', omega*rschw*r_calc / 1e5, '(km/s)'


    ! calculate the top values


    ! top flux from required accretion energy dissipation
    ! eddington accretion rate

    write (33,fmt33e)  'Eddington acc. rate', 24*3600.*365.*acc_edd / cgs_msun, '(Msun/yr)'
    write (33,fmt33e)  'Real acc. rate', acc_rate*24.*3600.*365.*acc_edd / cgs_msun, '(Msun/yr)'

    write (33,fmt33e)  'Energy flux', flux_acc, 'erg/cm2'
    write (33,fmt33f)  'Eff. temperature', temp_eff/1e6, 'MK'

    ! pressure from equation of state

    poly_k = ( (ro_0*0.01) / (miu*cgs_mhydr) * cgs_boltz * temp_eff     &
            + cgs_a / 3. * temp_eff ** 4 ) / (ro_0*0.01) ** gamm
    ! poly_k = 1e17


    ! disk height scale
    z_top_0 = z_top_start*zscale
    z_bottom = 100.


    ! allocate le memory
    allocate(z(ngrid), dz(ngrid), prs(ngrid), tau(ngrid), &
        alfven(ngrid),val(ngrid,n_vars),drv(ngrid,n_vars), csound1(ngrid), beta(ngrid), alfa_eff(ngrid))


    do iter = 1,iter_max

        print '(I6,A,Es12.4)', iter, ' poly_k', poly_k

        val(:,:) = 0.
        drv(:,:) = 0.
        val(1,v_ro) = ro_0
        val(1,v_flux) = 0.
        prs(1) = poly_k * ro_0**gamm
        val(1,v_pmag) =  prs(1) / beta_0
        val(1,v_temp) = findTemperature(prs(1),ro_0)
        tau(:) = 0.
        beta(1) = beta_0
        csound1(1) = sqrt(prs(1)/ro_0)
        last_i = ngrid

        ! disk height scale -- REPEATED
        zscale = sqrt( poly_k * ro_0 ** (gamm-1d0) ) / omega
        z_top_0 = z_top_start*zscale
        z_bottom = 100.
        ! generate the height scale
        z = (/ ( i / (ngrid-1.), i = 0,ngrid-1) /)
        z = exp( z*( log(z_top_0) - log(z_bottom) ) + log(z_bottom) - log(zscale) ) * zscale

        dz(2:ngrid) = z(2:ngrid) - z(1:ngrid-1)
        dz(1) = 2*dz(2) - dz(3)


        do i = 2, ngrid

            ! check if MRI takes place

            ! mri_active = 0.5*(1. + log(   (beta(i-1))/( 2 * csound1(i-1) / ( omega*r_calc*rschw ) )   )/log(1.1))
            ! mri_active = 0.5*(1. + log(   mri_limit(val(i-1,v_ro))/val(i-1,v_pmag)   )/log(1.+mri_tolerance))
            ! mri_active = real8_constraint(mri_active, min=0d0, max=1d0)
            mri_active = THR4POLY((log(val(i-1,v_pmag)) - log(val(i-1,v_ro)))/log(1.+mri_tolerance), 3)
            alfa_eff(i) = mri_active*alpha

            ! mri_active = l2i(mri_active_history .gt. 0)

            ! print '(I8,2Es12.4,2F12.4)', i, mri_limit(val(i-1,v_ro)), val(i-1,v_pmag), mri_active, mri_active_history

            ! RUNGE-KUTTA 4ord
            call compu_step(    z(i-1),         &
                                val(i-1,:),     &
                                (/ mri_active /),tmp_a)

            call compu_step(    z(i-1) + 0.5*dz(i),                 &
                                val(i-1,:) + 0.5*dz(i)*tmp_a,       &
                                (/ mri_active /), tmp_b)

            call compu_step(    z(i-1) + 0.5*dz(i),                 &
                                val(i-1,:) + 0.5*dz(i)*tmp_b,       &
                                (/ mri_active /), tmp_c)

            call compu_step(    z(i-1) + dz(i),                 &
                                val(i-1,:) + dz(i)*tmp_c,       &
                                (/ mri_active /), tmp_d)


            drv(i,:) = (tmp_a + 2*tmp_b + 2*tmp_c + tmp_d)/6.
            val(i,:) = drv(i,:) * dz(i) + val(i-1,:)

            ! ! if density goes negative we stop
            if ( val(i,v_ro) < 0 .or. ieee_is_nan(val(i,v_ro)) ) then
                last_i = i-1
                exit
            endif



            prs(i) = poly_k * val(i,v_ro)**gamm
            csound1(i) = sqrt( prs(i) / val(i,v_ro) )
            beta(i) = prs(i) / val(i,v_pmag)
            alfven(i) = sqrt( 2d0 * val(i,v_pmag) / val(i,v_ro) )


            tau(i) = tau(i-1) - val(i,v_ro) * cgs_kapes * dz(i)





        end do

        if ( iter .ne. iter_max ) then
            poly_k = poly_k * ( flux_acc / (val(last_i,v_flux) + 2. * val(last_i,v_pmag) * dzeta*omega*z(last_i) ) )**(2./3)
        endif
    end do

    !-------- opticale path

    tau = tau - tau(last_i)
    write (33,fmt33e)  'Deepest optical path: ', tau(1)
    ! write (33,'(A,2Es10.2)')  'Expected/computed beta0 ',  beta_0, prs(last_i)/prs_mag(last_i)
    !-------- save results

    write (33,fmt33i)  'Nr of iterations', iter

    open(file=outfn_ext, action='write', unit=10)
    write (10,fmt1s)  '0-Z', '1-TEMP', '2-RHO', '3-PRESS', '4-FLUX', '5-KAPPA', '6-TAU', &
        '7-FLUX_GEN', '8-BETA', '9-PMAG', '10-PB-MRI', '11-FLUXG_MAG', '12-PMAG_DERZ', '13-PRS_GAS',    &
        '14-PFLUX', '15-SOUNSP', '16-ALFV', '17-FF', '18-FLTOT', '19-Z/H', '20-tempful', &
         '21-tempflx', '22-temp-rad', '23-SS-LIM', '24-BET-LIM', '25-ALFAEFF', '26-FLEGN2'

    do i = 2, last_i, 4

        write (10,fmt1)  z(i),  prs(i)/val(i,v_ro) * miu * cgs_mhydr / cgs_boltz, val(i,v_ro), prs(i), val(i,v_flux),    &
            cgs_kapes, tau(i), drv(i,v_flux), beta(i), &
             val(i,v_pmag), mri_limit(val(i,v_ro)) ,   &
             alfa_eff(i)*omega*(prs(i)+val(i,v_pmag)) , drv(i,v_pmag), prs(i), 2*val(i,v_pmag)*z(i)*dzeta*omega, &
             csound1(i), alfven(i), omega*z(i), flux_acc, z(i)/zscale,     &
             findTemperature(prs(i),val(i,v_ro)), val(i,v_temp), (3. * prs(i) / cgs_a)**0.25,       &
             0.5 * beta(i) * omega * r_calc * rschw, 2 * csound1(i) / (omega * r_calc * rschw)

    end do


    write (33,fmt33e)  'Polytropic constant', poly_k, '[cgs]'
    write (33,fmt33f)  'Polytropic index', gamm
    write (33,fmt33f)  'Magnetic beta on equator', beta_0
    write (33,fmt33f)  'Magnetic alpha', alpha
    write (33,fmt33f)  'alpha_b / dzeta', alpha / dzeta
    write (33,fmt33f)  'Field rise dzeta parameter', dzeta
    write (33,fmt33e)  'Height scale',  zscale / 1e5, 'km'
    write (33,fmt33e)  'Top of calculations',  z_top_0 / 1e5, 'km'
    write (33,fmt33i)  'Nr of bins', ngrid
    write (33,fmt33e) 'dZ/Z', -log10(ngrid*1.) + log10(log(z_top_0) - log(z_bottom))

    close(unit=10)

    close(unit=33)


contains

    subroutine disk_read_local(cfg,errno)
        integer, intent(inout) :: errno
        type(confort_c), intent(inout) :: cfg

        call mincf_get(cfg, "z_top_start", buf, "12")
        read (buf,*) z_top_start

        call mincf_get(cfg,'beta_0',buf)
        read (buf,*) beta_0

        call mincf_get(cfg,'alpha',buf)
        read (buf,*) alpha

        call mincf_get(cfg,'rho_0',buf)
        read (buf,*) ro_0

        call MINCF_GET(cfg,'gamma',buf,'1.666666667')
        read (buf,*) gamm
    end subroutine

    elemental real(real64) function mri_limit(ro)
        real(real64), intent(in) :: ro
        mri_limit = 0.5 * sqrt( poly_k * ro**(gamm+1d0)        &
           * cgs_graw * (m_bh * cgs_msun) / (rschw*r_calc) )
    end function


    subroutine compu_step(z,val_in,params,der_out)
        real (dp), intent(in) :: z
        real (dp), intent(in) :: val_in(n_vars), params(n_params)
        real (dp), intent(out) :: der_out(n_vars)

        real(real64) :: pmag_der, mri_active, prs

        mri_active = params(p_mri)
        prs = poly_k * val_in(v_ro)**gamm

        pmag_der = (mri_active*(alpha / dzeta)*prs    &
                - (2.-mri_active*(alpha / dzeta))*val_in(v_pmag)) / z

        der_out(v_pmag) = pmag_der

        der_out(v_ro) = - val_in(v_ro)**(1-gamm) / (poly_k * gamm) * ( &
                    omega**2 * val_in(v_ro) * z + pmag_der  )

        der_out(v_flux) = -dzeta*omega*z*pmag_der

        der_out(v_temp) = - 3. * cgs_kapes * val_in(v_ro)    &
            / ( 16. * cgs_stef * (val_in(v_temp) ** 3) ) * val_in(v_flux)

    end subroutine


end program
