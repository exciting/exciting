!> Compute the Hilbert transform \(g(x)\) of a real function \(f(x)\) defined as
!> \[ g(x) = \frac{1}{\pi} {\rm p.v.} \int\limits_{x_{\rm min}}^{x_{\rm max}} \frac{f(x')}{x - x'} {\rm d}x' \]
!>
!> \(f(x)\) is expressed as piecewise cubic polynomials and the integrals are computed analytically for the cubic polynomials.
subroutine hilbert_transform( n, x, f, g )
  use precision, only: dp
  use constants, only: pi
  use asserts, only: assert
  !> number of sampling points
  integer, intent(in) :: n
  !> sampling points \(x_i\) in strictly ascending order
  real(dp), intent(in) :: x(n)
  !> input function \(f(x_i)\)
  real(dp), intent(in) :: f(n)
  !> Hilbert transform \(g(x_i)\)
  real(dp), intent(out) :: g(n)

  real(dp), parameter :: sixth = 1.0_dp / 6.0_dp

  integer :: i, j
  real(dp) :: d, dx, df

  real(dp), allocatable :: cf(:, :)

  ! get spline coefficients
  allocate( cf(3, n) )
  call spline( n, x, 1, f, cf )
  ! integrate
  g = 0.0_dp
  do i = 1, n
    do j = 1, n-1
      d = x(j+1) - x(j)
      dx = x(j) - x(i)
      df = f(j) - f(i)
      if (i == j) then
        g(i) = g(i) + sixth * d * (6*cf(1, j) + d * (3*cf(2, j) + d * 2*cf(3, j)))
      else if( i == j+1 ) then
        g(i) = g(i) + sixth * d * (6*cf(1, j) + d * (9*cf(2, j) + d * 11*cf(3, j)))
      else
        g(i) = g(i) + sixth * d * (6*cf(1, j) - 6*dx * (cf(2, j) - dx * cf(3, j)) + d * (3*cf(2, j) - 3*dx*cf(3, j) + 2*d*cf(3, j))) &
          + (df - dx * (cf(1, j) - dx * (cf(2, j) - dx * cf(3, j)))) * log( 1.0_dp + d / dx )
      end if
    end do
  end do
  g = g / pi

  deallocate( cf )
end subroutine hilbert_transform
