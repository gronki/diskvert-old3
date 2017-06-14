module query_param_binding

    use iso_c_binding
    use iso_fortran_env
    use globals
    use c_f_string_interop, only: c_f_string
    implicit none

contains

    function c_query_param_f(key_c) bind(C, name='dv_query_f') result(val)
        type(c_ptr), intent(in), value :: key_c
        real(c_float) :: val
        character(:), allocatable :: key

        call c_f_string(key_C, key)

        select case (key)
        case ("mass")
          val = m_bh
        case default
          error stop "no such key"
        end select

    end function

end module
