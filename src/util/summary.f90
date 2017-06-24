module diskvert_summary

    character (len=*) , parameter :: fmt33e = "(A30,Es15.4,1X,A8)"
    character (len=*) , parameter :: fmt33s = "(A30,A15,1X,A8)"
    character (len=*) , parameter :: fmt33f = "(A30,F15.4,1X,A8)"
    character (len=*) , parameter :: fmt33i = "(A30,I15,1X,A8)"
    character (len=*) , parameter :: fmtdescr = '(A18,I4,A12,I10,I6)'

    integer :: length_label = 20
    integer :: length_value = 12
    integer :: length_unit = 10
    integer :: length_fltprec = 4


    interface SUMMPRINT
        module procedure SUMMPRINT_I
        module procedure SUMMPRINT_S
        module procedure SUMMPRINT_F
    end interface



contains

    subroutine assertion_old(l,m,die)
        logical, intent(in) :: l
        logical, intent(in),optional :: die
        character(len=*), intent(in) :: m
        if(.not.(l)) write (0,*) m
    end subroutine


    subroutine SUMMPRINT_S(file,label,strval,unit)
        character(len=*), intent(in) :: label
        character(len=*), intent(in) :: strval
        character(len=*), intent(in), optional :: unit
        integer, intent(in), optional :: file
        character(len=128) :: buf
        character(len=128) :: fmt
        integer :: ios
        integer :: funit

        funit = 6
        if ( PRESENT(file) )    funit=file

        if ( PRESENT(unit) ) then
            write(buf,'(A)') ADJUSTL(unit)
        else
            write(buf,'(A)')  ' '
        end if

        call summformat_s(fmt)
        write(unit=funit, fmt=fmt, iostat=ios) label,strval,buf(1:length_unit)

        call assertion_old(ios.eq.0, 'could not write to given unit', DIE=.false.)

    end subroutine

    subroutine SUMMPRINT_I(file,label,intval,unit)
        character(len=*), intent(in) :: label
        integer, intent(in) :: intval
        character(len=*), intent(in), optional :: unit
        integer, intent(in), optional :: file
        character(len=128) :: buf
        character(len=128) :: fmt
        integer :: ios
        integer :: funit

        funit = 6
        if ( PRESENT(file) )    funit=file

        if ( PRESENT(unit) ) then
            write(buf,'(A)') ADJUSTL(unit)
        else
            write(buf,'(A)')  ' '
        end if

        call summformat_i(fmt)
        write(unit=funit, fmt=fmt, iostat=ios) label,intval,buf(1:length_unit)

        call assertion_old(ios.eq.0, 'could not write to given unit', DIE=.false.)

    end subroutine


    subroutine SUMMPRINT_F(file,label,fltval,unit)
        character(len=*), intent(in) :: label
        real, intent(in) :: fltval
        character(len=*), intent(in), optional :: unit
        integer, intent(in), optional :: file
        character(len=128) :: buf, buf2
        character(len=128) :: fmt
        integer :: ios
        integer :: funit

        funit = 6
        if ( PRESENT(file) )    funit=file

        if ( PRESENT(unit) ) then
            write(buf,'(A)') ADJUSTL(unit)
        else
            write(buf,'(A)')  ' '
        end if

        call summformat_f(fmt,fltval)
        write(unit=funit, fmt=fmt, iostat=ios) label,fltval,buf(1:length_unit)

        call assertion_old(ios.eq.0, 'could not write to given unit', DIE=.false.)

    end subroutine

    !!!!!!!!!!!!!!! SUMMFORMAT !!!!!!!!!!!!!!!

    subroutine summformat_f(fmt, fltval)
        character(len=*), intent(out) :: fmt
        real, intent(in), optional :: fltval
        character(len=32) :: valfmt
        real :: fv
        if ( present(fltval) ) then
            fv = fltval
        else
            fv = 1e7
        end if

        if ( abs(fv) .gt. 0.01 .and. abs(fv) .lt. 100. ) then
            write(valfmt,'(A,I0,A,I0)') 'F', length_value, '.', length_fltprec
        else
            write(valfmt,'(A,I0,A,I0)') 'Es', length_value, '.', length_fltprec
        end if

        write(fmt,'(A,I0,A,A,A,I0,A)') '(A', length_label, ',', valfmt,  ',A',length_unit+1,')'
    end subroutine

    subroutine summformat_i(fmt)
        character(len=*), intent(out) :: fmt
        character(len=32) :: valfmt
        write(valfmt,'(A,I0)') 'I', length_value
        write(fmt,'(A,I0,A,A,A,I0,A)') '(A', length_label, ',', valfmt,  ',A',length_unit+1,')'
    end subroutine

    subroutine summformat_s(fmt)
        character(len=*), intent(out) :: fmt
        character(len=32) :: valfmt
        write(valfmt,'(A,I0)') 'A', length_value
        write(fmt,'(A,I0,A,A,A,I0,A)') '(A', length_label, ',', valfmt,  ',A',length_unit+1,')'
    end subroutine


end module
