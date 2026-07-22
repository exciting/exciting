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
    !> currently implemented: 'lorentzian' (default), 'gaussian'
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
    real(kind=dp) :: w
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
        w = width
        p0 => poly_lorentz_0
        select case (mde)
          case ('linear')
            pn => poly_lorentz_1
          case ('cubic')
            pn => poly_lorentz_3
        end select
      case ('gaussian')
        w = width / sqrt(log(2.0_dp))
        p0 => poly_gauss_0
        select case (mde)
          case ('linear')
            pn => poly_gauss_1
          case ('cubic')
            pn => poly_gauss_3
        end select
      case default
        error stop '('//PROCEDURE_NAME//') ' // &
          'Invalid convolution kernel `'//ker//'`.'
    end select

    do i = 1, n
      y(i) = y(1)*p0( c(:, 1), w, x(i), x(1), x(1) ) + y(n)*(1.0_dp - p0( c(:, n), w, x(i), x(n), x(n) ))
      do j = 1, n-1
        y(i) = y(i) + (pn( c(:, j), w, x(i), x(j), x(j+1) ) - pn( c(:, j), w, x(i), x(j), x(j) ))
      end do
    end do

  contains

    ! ****************************
    ! Lorentzian kernel
    !
    ! \int_{-\infty}^x_0 K(t-x,w)\, {\rm d}t
    function poly_lorentz_0( c, w, x, x0, t ) result( p )
      real(kind=dp), intent(in) :: c(0:), w, x, x0, t
      real(kind=dp) :: p

      p = (0.5_dp + atan( (x0 - x) / w ))
      p = p / pi
    end function

    ! \int (c_0 + c_1(t-x_0))\, K(t-x,w)\, {\rm d}t
    function poly_lorentz_1( c, w, x, x0, t ) result( p )
      real(kind=dp), intent(in) :: c(0:), w, x, x0, t
      real(kind=dp) :: p

      p = 2*(c(0) + c(1)*(x - x0))*atan( (t - x) / w ) + &
          w*c(1)*log( w**2 + (t - x)**2 )
      p = p / (2*pi)
    end function

    ! \int (c_0 + c_1(t-x_0) + c_2(t-x_0)^2 + c_3(t-x_0)^3)\, K(t-x,w)\, {\rm d}t
    function poly_lorentz_3( c, w, x, x0, t ) result( p )
      real(kind=dp), intent(in) :: c(0:), w, x, x0, t
      real(kind=dp) :: p

      p = 2*w*(t - x)*(c(2) + 3*c(3)*(x - x0)) + &
          w*c(3)*(t - x)**2 + &
          2*(c(0) - w**2*c(2) + (x - x0)*(c(1) - 3*w**2*c(3) + (x - x0)*(c(2) + c(3)*(x - x0))))*atan( (t - x) / w ) + &
          w*(c(1) + 2*c(2)*(x - x0) - c(3)*(w**2 - 3*(x - x0)**2))*log( w**2 + (t - x)**2 )
      p = p / (2*pi)
    end function
    ! ****************************

    ! ****************************
    ! Gaussian kernel
    !
    function poly_gauss_0( c, w, x, x0, t ) result( p )
      real(kind=dp), intent(in) :: c(0:), w, x, x0, t
      real(kind=dp) :: p

      p = 1.0_dp + erf( (x0 - x) / w )
      p = p / 2
    end function

    function poly_gauss_1( c, w, x, x0, t ) result( p )
      real(kind=dp), intent(in) :: c(0:), w, x, x0, t
      real(kind=dp) :: p

      p = (c(0) + c(1)*(x - x0))*erf( (t - x) / w ) - &
          c(1)*w/sqrt(pi)*exp( -( (t - x) / w )**2 )
      p = p / 2
    end function

    function poly_gauss_3( c, w, x, x0, t ) result( p )
      real(kind=dp), intent(in) :: c(0:), w, x, x0, t
      real(kind=dp) :: p

      p = (2*c(0) + w**2*c(2) + (x - x0)*(2*c(1) + 3*w**2*c(3) + 2*(x - x0)*(c(2) + c(3)*(x - x0))))*erf( (t - x) / w ) - &
          2*(c(1) + c(2)*(t + x - 2*x0) + c(3)*(w**2 + t**2 + t*x + x**2 - 3*(t + x - x0)*x0))*w/sqrt(pi)*exp( -( (t - x) / w )**2 ) 
      p = p / 4
    end function
    ! ****************************
  end subroutine smoothen

end module convolution
