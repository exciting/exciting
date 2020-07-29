    
subroutine get_selfc(n, x, y, x0, f, df)
    use mod_pade
    implicit none
    integer,    intent(in)  :: n
    real(8),    intent(in)  :: x(n)
    complex(8), intent(in)  :: y(n)
    real(8),    intent(in)  :: x0
    complex(8), intent(out) :: f
    complex(8), intent(out) :: df
    ! local
    integer(4) :: i, i0, ii, j
    integer(4) :: np, np2
    real(8) :: xx, t1, t2
    real(8),    allocatable :: xa(:), c(:)
    complex(8), allocatable :: ya(:)
    real(8),    external    :: polynom

    ! only 1 point
    if (n <= 4) then
        write(*,*) 'ERROR(get_selfc): Not enough data points!'
        stop
    end if

    xx = x0
    if ( xx < x(1) ) xx = x(1)
    if ( xx > x(n) ) xx = x(n)

    !-----------------------
    ! polynomial fitting
    !-----------------------
    np = 4 ! cubic polynomial fitting
    np2 = np/2
    allocate(xa(np), ya(np), c(np))
    do i = 1, n
        if (x(i) >= xx) then
            if (i <= np2) then
                i0 = 1
            else if (i > n-np2) then
                i0 = n-np+1
            else
                i0 = i-np2
            end if
            do j = 1, np
                ii = i0+j-1
                xa(j) = x(ii)
                ya(j) = y(ii)
            end do
            t1 = polynom(0, np, xa, dble(ya), c, xx)
            t2 = polynom(0, np, xa, aimag(ya), c, xx)
            f  = cmplx(t1, t2, 8)
            t1 = polynom(1, np, xa, dble(ya), c, xx)
            t2 = polynom(1, np, xa, aimag(ya), c, xx)
            df = cmplx(t1, t2, 8)
            ! call pade_approximant(np, cmplx(xa,0.d0,8), ya, cmplx(xx,0.d0,8), f, df)
            exit
        end if
    end do
    deallocate(xa, ya, c)

    return
end subroutine
