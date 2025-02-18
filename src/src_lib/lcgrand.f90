!> Simple compiler independent random number generator.
!> Generate pseudo random numbers in the interval [0,1) using a 
!> [linear congruential generator](https://en.wikipedia.org/wiki/Linear_congruential_generator).
!> 
!> Source: Numerical Recipes from the **quick and dirty generators** list, Chapter 7.1, Eq. 7.1.6
!> parameters from Knuth and H. W. Lewis
!
! MRM (2025): Using long integers. Cray fails otherwise.
!
subroutine lcgrand( v, n, s)
  use precision, only: dp, long_int
  !> vector to be filled with random numbers
  real(dp), intent(out) :: v(*)
  !> number or random numbers
  integer, intent(in) :: n
  !> seed for the generator
  integer, intent(in) :: s

  integer(long_int), parameter :: m = 2_long_int**32_long_int       !! modulus
  integer(long_int), parameter :: a = 1664525_long_int     !! multiplier
  integer(long_int), parameter :: c = 1013904223_long_int  !! increment

  integer(dp) :: i, j, k

  j = s
  do i = 1, n
    k = mod( a*j + c, m)
    v(i) = real(k, kind=dp)/real(m, kind=dp)
    j = k
  end do
end subroutine
