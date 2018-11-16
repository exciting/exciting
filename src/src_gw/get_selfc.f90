    
subroutine get_selfc(n, x, y, x0, f, df)
    implicit none
    integer,    intent(in)  :: n
    real(8),    intent(in)  :: x(n)
    complex(8), intent(in)  :: y(n)
    real(8),    intent(in)  :: x0
    complex(8), intent(out) :: f
    complex(8), intent(out) :: df
    ! local
    integer    :: i, i0, ii, j
    integer    :: np, np2
    real(8)    :: dx, t1, t2, xx
    real(8),    allocatable :: c(:)
    complex(8), allocatable :: ya(:)
    real(8),    external    :: polynom

    ! only 1 point
    if (n <= 1) then
        write(*,*) 'ERROR(get_selfc): Not enough data points!'
        stop
    end if

    xx = max(x0, x(1))
    xx = min(x0, x(n))

    !-----------------------
    ! polynomial fitting
    !-----------------------
    np = 8
    np2 = np/2
    allocate(ya(np), c(np))
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
                ya(j) = y(ii)
            end do
            t1 = polynom( 0, np, x(i0), dble(ya), c, xx)
            t2 = polynom( 0, np, x(i0), aimag(ya), c, xx)
            f  = cmplx(t1, t2, 8)
            t1 = polynom( 1, np, x(i0), dble(ya), c, xx)
            t2 = polynom( 1, np, x(i0), aimag(ya), c, xx)
            df = cmplx(t1, t2, 8)
            exit
        end if
    end do
    deallocate(ya, c)

    return
end subroutine
