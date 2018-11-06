
subroutine clenshaw_curtis_compute ( order, x, w )

!*****************************************************************************80
!
!! CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis quadrature rule.
!
!  Discussion:
!
!    Our convention is that the abscissas are numbered from left to right.
!
!    The rule is defined on [-1,1].
!
!    The integral to approximate:
!
!      Integral ( -1 <= X <= 1 ) F(X) dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!    1 <= ORDER.
!
!    Output, real ( kind = 8 ) X(ORDER), the abscissas.
!
!    Output, real ( kind = 8 ) W(ORDER), the weights.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) theta
  real ( kind = 8 ) w(order)
  real ( kind = 8 ) x(order)

  if ( order < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CLENSHAW_CURTIS_COMPUTE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', order
    stop
  end if

  if ( order == 1 ) then
    x(1) = 0.0D+00
    w(1) = 2.0D+00
    return
  end if

  do i = 1, order
    x(i) = cos ( real ( order - i, kind = 8 ) * r8_pi &
               / real ( order - 1, kind = 8 ) )
  end do

  x(1) = -1.0D+00
  if ( mod ( order, 2 ) == 1 ) then
    x((order+1)/2) = 0.0D+00
  end if
  x(order) = +1.0D+00

  do i = 1, order

    theta = real ( i - 1, kind = 8 ) * r8_pi &
          / real ( order - 1, kind = 8 )

    w(i) = 1.0D+00

    do j = 1, ( order - 1 ) / 2

      if ( 2 * j == ( order - 1 ) ) then
        b = 1.0D+00
      else
        b = 2.0D+00
      end if

      w(i) = w(i) - b * cos ( 2.0D+00 * real ( j, kind = 8 ) * theta ) &
           / real ( 4 * j * j - 1, kind = 8 )

    end do

  end do

  w(1)         =           w(1)         / real ( order - 1, kind = 8 )
  w(2:order-1) = 2.0D+00 * w(2:order-1) / real ( order - 1, kind = 8 )
  w(order)     =           w(order)     / real ( order - 1, kind = 8 )

  return
end

subroutine rescale(a, b, n, x, w)

!*****************************************************************************80
!
!! RESCALE rescales a quadrature rule from [-1,+1] to [A,B].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the new interval.
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input/output, real ( kind = 8 ) X(N), on input, the abscissas for [-1,+1].
!    On output, the abscissas for [A,B].
!
!    Input/output, real ( kind = 8 ) W(N), on input, the weights for [-1,+1].
!    On output, the weights for [A,B].
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  x(1:n) = ( ( a + b ) + ( b - a ) * x(1:n) ) / 2.0D+00
  w(1:n) = ( b - a ) * w(1:n) / 2.0D+00

  return
end
