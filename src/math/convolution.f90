!> Convolutions of 1D functions.
module convolution
  use precision, only: dp
  implicit none
  private

  public :: smoothen
  
contains

  !> Smoothen function by convolution with selected kernel of given width.
  !> The width is given as half width at half maximum (HWHM).
  subroutine smoothen( x, y, width, kernel, mode )
    use constants, only: pi
    !> x-values functions is given at (in increasing order)
    real(kind=dp), intent(in) :: x(:)
    !> on input: function values at given x   
    !> on output: convoluted function
    real(kind=dp), intent(inout) :: y(:)
    !> kernel width
    real(kind=dp), intent(in) :: width
    !> convolution kernel   
    !> currently implemented: 'lorentzian' (default)
    character(len=*), optional, intent(in) :: kernel
    !> integration mode   
    !> currently implemented: 'linear' (default), 'cubic'
    character(len=*), optional, intent(in) :: mode

    abstract interface
      function poly_integral( c, w, x, x0, t ) result( p )
        use precision, only: dp
        real(kind=dp), intent(in) :: c(0:), w, x, x0, t
        real(kind=dp) :: p
      end function
    end interface

    character(len=*), parameter :: PROCEDURE_NAME = 'smoothen'
  
    character(len=:), allocatable :: ker, mde
    integer :: n, i, j
    real(kind=dp), allocatable :: c(:,:)
    procedure(poly_integral), pointer :: p0, pn

    ker = 'lorentzian'
    if (present(kernel)) ker = trim( adjustl( kernel ) )
    mde = 'linear'
    if (present(mode)) mde = trim( adjustl( mode ) )

    n = size( x, dim=1 )
    if (size( y, dim=1 ) /= n) error stop '('//PROCEDURE_NAME//') ' // &
      'Arrays `x` and `y` must be of same length.'
    if (width < epsilon(1.0_dp)) error stop '('//PROCEDURE_NAME//') ' // &
      '`width` must be positive.'

    ! set polynomial coefficients
    allocate( c(0:3, n), source=0.0_dp )
    c(0, :) = y
    select case (mde)
      case ('linear')
        do i = 1, n-1
          c(1, i) = (y(i+1) - y(i)) / (x(i+1) - x(i))
        end do
      case ('cubic')
        call spline( n, x, 1, y, c(1:3, :) )
      case default
        error stop '('//PROCEDURE_NAME//') ' // &
          'Invalid integration mode `'//mde//'`.'
    end select

    ! set polynomial integrals
    select case (ker)
      case ('lorentzian')
        p0 => poly_lorentz_0
        select case (mde)
          case ('linear')
            pn => poly_lorentz_1
          case ('cubic')
            pn => poly_lorentz_3
        end select
      case default
        error stop '('//PROCEDURE_NAME//') ' // &
          'Invalid convolution kernel `'//ker//'`.'
    end select

    do i = 1, n
      y(i) = 0.5_dp*(y(1) + y(n)) + p0( c(:, 1), width, x(i), x(1), x(1) ) - p0( c(:, n), width, x(i), x(n), x(n) )
      do j = 1, n-1
        y(i) = y(i) + (pn( c(:, j), width, x(i), x(j), x(j+1) ) - pn( c(:, j), width, x(i), x(j), x(j) ))
      end do
    end do

  contains

    ! ****************************
    ! Lorentzian kernel
    !
    function poly_lorentz_0( c, w, x, x0, t ) result( p )
      real(kind=dp), intent(in) :: c(0:), w, x, x0, t
      real(kind=dp) :: p

      p = c(0)*atan( (t - x) / w )
      p = p / pi
    end function

    function poly_lorentz_1( c, w, x, x0, t ) result( p )
      real(kind=dp), intent(in) :: c(0:), w, x, x0, t
      real(kind=dp) :: p

      p = 4*(c(0) + c(1)*(x - x0))*atan( (t - x) / w ) + &
          2*w*c(1)*log( w**2 + (t - x)**2 )
      p = p / (4*pi)
    end function

    function poly_lorentz_3( c, w, x, x0, t ) result( p )
      real(kind=dp), intent(in) :: c(0:), w, x, x0, t
      real(kind=dp) :: p

      p = 16*w*(t - x0)*(c(2) + 2*c(3)*(x - x0)) + &
          8*w*c(3)*(t - x0)**2 + &
          16*(c(0) - w**2*c(2) + (x - x0)*(c(1) - 3*w**2*c(3) + (x - x0)*(c(2) + c(3)*(x - x0))))*atan( (t - x) / w ) + &
          8*w*(c(1) + 2*c(2)*(x - x0) - c(3)*(w**2 - 3*(x - x0)**2))*log( w**2 + (t - x)**2 )
      p = p / (16*pi)
    end function
    ! ****************************
  end subroutine smoothen

end module convolution
