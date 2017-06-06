module slf_space

    use ieee_arithmetic
    use iso_fortran_env

    use precision

    implicit none

    interface space
        module procedure space1, space2, space3
    end interface

contains


    pure subroutine space1(x)

        real(fp), dimension(:), intent(inout) :: x
        integer :: i,n

        n = size(x)

        forall (i = 1:n)
            x(i) = real(i - 1, fp) / ( n - 1 )
        end forall

    end subroutine


    pure subroutine space2(x,dim)
        integer, intent(in) :: dim
        real(fp), dimension(:,:), intent(inout) :: x
        integer :: i,n

        if ( dim == 1 ) then
            n = size(x,1)
            forall (i = 1:n)
                x(i,:) = real(i - 1, fp) / ( n - 1 )
            end forall
        elseif ( dim == 2 ) then
            n = size(x,2)
            forall (i = 1:n)
                x(:,i) = real(i - 1, fp) / ( n - 1 )
            end forall
        end if

    end subroutine


    pure subroutine space3(x,dim)
        integer, intent(in) :: dim
        real(fp), dimension(:,:,:), intent(inout) :: x
        integer :: i,n

        if ( dim == 1 ) then
            n = size(x,1)
            forall (i = 1:n)
                x(i,:,:) = real(i - 1, fp) / ( n - 1 )
            end forall
        elseif ( dim == 2 ) then
            n = size(x,2)
            forall (i = 1:n)
                x(:,i,:) = real(i - 1, fp) / ( n - 1 )
            end forall
        elseif ( dim == 3 ) then
            n = size(x,3)
            forall (i = 1:n)
                x(:,:,i) = real(i - 1, fp) / ( n - 1 )
            end forall
        end if

    end subroutine


end module
