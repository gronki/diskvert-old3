module c_f_string_interop

    use iso_c_binding

contains

    ! source: https://stackoverflow.com/a/11443635

    subroutine c_f_string(c_string, f_string)
      use, intrinsic :: iso_c_binding, only: c_char, c_null_char, c_f_pointer
      !----
      type(c_ptr), intent(in) :: c_string
      character(:), intent(out), allocatable :: f_string
      !----
      ! Array for accessing string pointed at by C pointer
      character(kind=c_char), pointer :: string_ptr(:)
      integer :: i    ! string index
      interface
        ! Steal std C library function rather than writing our own.
        function strlen(s) bind(c, name='strlen')
          use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
          implicit none
          !----
          type(c_ptr), intent(in), value :: s
          integer(c_size_t) :: strlen
        end function strlen
      end interface
      !****
      ! Map C pointer to fortran character array
      call c_f_pointer(c_string, string_ptr, [strlen(c_string)])
      ! Allocate fortran character variable to the c string's length
      allocate(character(size(string_ptr)) :: f_string)
      ! Copy across (with possible kind conversion) characters
      forall (i = 1:size(string_ptr)) f_string(i:i) = string_ptr(i)
    end subroutine c_f_string   

end module
