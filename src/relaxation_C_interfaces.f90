module relaxation_C_interfaces

    use iso_c_binding
    use iso_fortran_env, only: r64 => real64
    
    use relaxation

contains


    subroutine generate_coefficients_c_interface(X,nx,Y,ny,A,na,M) &
            & bind(C, name='generate_coefficients')

        integer(c_int), intent(in), value :: nx,ny,na
        real(c_double), intent(in), dimension(nx) :: X
        real(c_double), intent(in), dimension(nx*ny) :: Y
        real(c_double), intent(out), dimension(nx*na) :: A
        real(c_double), intent(out), dimension(nx*na,nx*ny) :: M

        type(model_t) :: model

        call model % init(compton = .false., magnetic = .false., conduction = .false.)
        call model % matrix(x,Y,M,A)

    end subroutine

    subroutine generate_corrections_c_interface(X,nx,Y,dY,ny) &
            & bind(C, name='generate_corrections')

        integer(c_int), intent(in), value :: nx,ny
        real(c_double), intent(in), dimension(nx) :: X
        real(c_double), intent(in), dimension(nx*ny) :: Y
        real(c_double), intent(out), dimension(nx*ny) :: dY
        real(r64) :: M(nx*ny,nx*ny)

        type(model_t) :: model

        call model % init(compton = .false., magnetic = .false., conduction = .false.)
        call model % advance(x,Y,M,dY)

    end subroutine


end module
