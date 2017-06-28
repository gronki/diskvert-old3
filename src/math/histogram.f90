module slf_histo

    use iso_fortran_env, only: r64 => real64
    implicit none

contains

    subroutine histo(x,min,max,h,loc)
        real(r64), intent(in) :: x(:)
        real(r64), intent(in) :: min,max
        real(r64), intent(out) :: h(:)
        real(r64), intent(out), optional :: loc(:)
        real(r64) :: bin
        integer :: n,nb,i,k
        n = size(x)
        nb = size(h)
        bin = (max-min)/nb
        h = 0
        do i=1,n
            k = floor((x(i) - min) / bin) + 1
            if (k .ge. 1 .and. k .le. nb)  h(k) = h(k) + 1d0 / n
        end do
        if (present(loc)) then
            do i=1,nb
                loc(i) = (i-0.5) * bin + min
            end do
        end if
    end subroutine

    subroutine histo2d(x,y,xlim,ylim,h,locx,locy)
        real(r64), intent(in) :: x(:), y(:)
        real(r64), intent(in) :: xlim(2), ylim(2)
        real(r64), intent(out) :: h(:,:)
        real(r64), intent(out), optional :: locx(:,:), locy(:,:)
        real(r64) :: xbin, ybin
        integer :: n,nx,ny,i,j,k
        n = size(x)
        nx = size(h,1)
        ny = size(h,2)
        xbin = (xlim(2)-xlim(1))/nx
        ybin = (ylim(2)-ylim(1))/ny
        h = 0
        do i=1,n
            j = floor((x(i) - xlim(1)) / xbin) + 1
            k = floor((y(i) - ylim(1)) / ybin) + 1
            if (j .ge. 1 .and. j .le. nx        &
                .and. k .ge. 1 .and. k .le. ny)  h(j,k) = h(j,k) + 1d0 / n
        end do
        if (present(locx)) then
            do i=1,nx
                do j=1,ny
                    locx(i,j) = (i-0.5) * xbin + xlim(1)
                end do
            end do
        end if
        if (present(locy)) then
            do i=1,nx
                do j=1,ny
                    locy(i,j) = (j-0.5) * ybin + ylim(1)
                end do
            end do
        end if
    end subroutine


end module
