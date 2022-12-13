!> Simple compiler independent random number generator.
!> Generate pseudo random numbers in the interval [0,1) using a 
!> [linear congruential generator](https://en.wikipedia.org/wiki/Linear_congruential_generator).
!> 
!> Source: Numerical Recipes from the **quick and dirty generators** list, Chapter 7.1, Eq. 7.1.6
!> parameters from Knuth and H. W. Lewis
subroutine lcgrand( v, n, s)
  use precision, only: dp
  !> vector to be filled with random numbers
  real(dp), intent(out) :: v(*)
  !> number or random numbers
  integer, intent(in) :: n
  !> seed for the generator
  integer, intent(in) :: s

  integer(dp), parameter :: m = 2**32       !! modulus
  integer(dp), parameter :: a = 1664525     !! multiplier
  integer(dp), parameter :: c = 1013904223  !! increment

  integer(dp) :: i, j, k

  j = s
  do i = 1, n
    k = mod( a*j + c, m)
    v(i) = dble(k)/dble(m)
    j = k
  end do
end subroutine
