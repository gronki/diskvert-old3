module slf_random

    use slf_interpol
    use iso_fortran_env

    implicit none

    character(len=*), parameter, private :: debugfmt = '(A,":",I0,": ",A)'

contains

    subroutine random_init_urandom()
        integer :: n,err
        integer, allocatable :: seed(:)

        call random_seed(size=n)
        allocate(seed(n))
        if (.not. allocated(seed)) then
            write (0,*) 'init_random: allocation failed'
            return
        end if
        open(unit=33,file='/dev/urandom',action='read',form='unformatted',iostat=err)
        if ( err .ne. 0 ) then
            write(0,*) 'init_random: /dev/urandom not found'
            deallocate(seed)
            return
        end if
        read (33) seed
        close(33)
        call random_seed(put=seed)
        deallocate(seed)
    end subroutine

    subroutine random_pdf(x,pdf,xout,cdf)
        real(real64), intent(in) :: x(:)
        real(real64), intent(out) :: xout(:)
        real(real64), intent(in) :: pdf(size(x))
        real(real64), intent(out) :: cdf(size(x))
        real(real64) :: rand(size(xout)), dx(size(x))
        integer :: err,i

        if(size(x) .le. 3) write (error_unit,debugfmt) __FILE__, __LINE__, &
                & 'rand_pdf: given PDF has too few points'

        do i=1,size(x)
            if ( i .gt. 1 ) then
                dx(i) = x(i) - x(i-1)
            else
                dx(i) = x(2) - x(1)
            end if
        end do

        if (.not.all(pdf .ge. 0)) write (error_unit,debugfmt) __FILE__, __LINE__, &
                & 'Found PDF<0'

        cdf(1) = pdf(1)*dx(1)
        if ( cdf(1) .lt. epsilon(1d0) ) cdf(1) = epsilon(1d0)

        do i = 2,size(x)
            cdf(i) = cdf(i-1) + pdf(i)*dx(i)
            if ( cdf(i) .le. cdf(i-1) ) then
                cdf(i) = cdf(i-1) * ( 1 + epsilon(1d0) )
            end if
        end do

        cdf = cdf / cdf(size(x))

        call random_number(rand)

        do i=1,size(rand)
            call interpol(cdf,x,rand(i),xout(i))
        end do

    end subroutine

end module
