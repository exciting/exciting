!> Savitzky-Golay filter for smoothing and derivatives of real 1D functions.
module savitzky_golay
  use precision, only: dp
  implicit none
  private

  character(len=*), private, parameter :: MODULE_NAME = 'savitzky_golay'

  public :: savgol

  !> Savitzky-Golay filter
  interface savgol
    module procedure :: savgol_many, savgol_single
  end interface savgol
  
contains

  !> Computes a smooth approximation to the function \(y = f(x)\) given by
  !> noisy data points \((x_i,y_i)\) by locally fitting a polynomial of given
  !> degree to the data points \(x_{i-m},\ldots,x_i,\ldots,x_{i+m}\) in a 
  !> least squares fashion and evaluating the polynomial at \(x_i\).
  subroutine savgol_many( x, y, degree, window, derivative, weighting )
    use grid_utils, only: linspace
    !> values \(x_i\) in increasing order
    real(kind=dp), intent(in) :: x(:)
    !> on input: values \(y_i\) in first column   
    !> on output: smooth approximation \(f(x_i)\) or selected derivatives \(f^{(d)}(x_i)\)
    real(kind=dp), intent(inout) :: y(:,:)
    !> polynomial degree 
    !>  (default: `3`)
    integer, optional, intent(in) :: degree
    !> window length \(m\) (\(2m+1\) points centered around \(x_i\) are considered for least squares fit) 
    !>  (default: `3`)
    integer, optional, intent(in) :: window
    !> order of derivatives to compute (derivatives of order greater than polynomial degree are zero) 
    !>  (default: `[0]` = smooth function only )
    integer, optional, intent(in) :: derivative(:)
    !> weighting function   
    !>  valid options: `'rectangle'` (no weighting)   
    !>                 `'hann'` (Hann function, default)
    character(len=*), optional, intent(in) :: weighting 
  
    character(len=*), parameter :: PROCEDURE_NAME = 'savgol_many'

    integer :: n, deg, win, d
    logical :: equidist
    character(len=:), allocatable :: wfun

    integer, allocatable :: dord(:)
    real(kind=dp), allocatable :: x0(:,:), s(:)

    ! SET DEFAULTS
    ! polynomial degree
    deg = 3
    if (present(degree)) deg = degree
    ! window length
    win = 3
    if (present(window)) win = window
    ! derivatives
    dord = [0]
    if (present(derivative)) dord = derivative
    ! weighting function
    wfun = 'hann'
    if (present(weighting)) wfun = weighting
    ! number of data points
    n = size(x)

    ! CHECK INPUT
    if (deg < 0) error stop 'Error ('//MODULE_NAME//':'//PROCEDURE_NAME//'): '// &
      'Polynomial degree `deg` must not be negative.'
    if (2*win < deg) error stop 'Error ('//MODULE_NAME//':'//PROCEDURE_NAME//'): '// &
      'Window length too small for given polynomial degree. It must hold `2*win >= deg`.'
    if (n < 2*win+1) error stop 'Error ('//MODULE_NAME//':'//PROCEDURE_NAME//'): '// &
      'Number of data points `n` too small for given window length `win`. It must hold `n >= 2*win+1`.'
    if (size(y, dim=1) < n) error stop 'Error ('//MODULE_NAME//':'//PROCEDURE_NAME//'): '// &
      'First dimension of `y` is too small.'
    if (size(y, dim=2) < size(dord)) error stop 'Error ('//MODULE_NAME//':'//PROCEDURE_NAME//'): '// &
      'Second dimension of `y` is too small.'
    if (any(dord < 0)) error stop 'Error ('//MODULE_NAME//':'//PROCEDURE_NAME//'): '// &
      'Only derivatives of order >= 0 allowed.'
    if (wfun /= 'rectangle' .and. wfun /= 'hann') error stop 'Error ('//MODULE_NAME//':'//PROCEDURE_NAME//'): '// &
      'Invalid weighting function.'
    if (any( x(2:n) - x(1:n-1) < 2*epsilon(x)*abs(x(n) - x(1)) )) error stop 'Error ('//MODULE_NAME//':'//PROCEDURE_NAME//'): '// &
      '`x` values are not strictly increasing.'

    ! check if sampling is equidistant
    equidist = all( abs( x(2:n) - x(1:n-1) - (x(2) - x(1)) ) < 2*epsilon(x)*(x(n) - x(1)) )
    if (.not. equidist .and. any(dord > 1)) error stop 'Error ('//MODULE_NAME//':'//PROCEDURE_NAME//'): '// &
      'Only 0th and 1st derivative supported for non-equidistant sampling.'

    if (equidist) then
      call savgol_equi( x, y, deg, win, dord, wfun )
    else
      s = linspace( x(1), x(n), n )
      call savgol_equi( s, y, deg, win, dord, wfun )
      if (any(dord == 1)) then
        allocate( x0(n, 1) )
        x0(:, 1) = x
        call savgol_equi( s, x0, deg, win, [1], wfun )
        do d = 1, size(dord)
          if (dord(d) /= 1) cycle
          y(:, d) = y(:, d) / x0(:, 1) 
        end do
      end if
    end if
  end subroutine savgol_many

  !> See [[savgol_many(subroutine)]].
  subroutine savgol_single( x, y, degree, window, derivative, weighting )
    !> values \(x_i\) in increasing order
    real(kind=dp), intent(in) :: x(:)
    !> on input: values \(y_i\)  
    !> on output: smooth approximation \(f(x_i)\) or selected derivative \(f^{(d)}(x_i)\)
    real(kind=dp), target, intent(inout) :: y(:)
    !> polynomial degree 
    !>  (default: `3`)
    integer, optional, intent(in) :: degree
    !> window length \(m\) (\(2m+1\) points centered around \(x_i\) are considered for least squares fit) 
    !>  (default: `3`)
    integer, optional, intent(in) :: window
    !> order of derivative to compute (derivatives of order greater than polynomial degree are zero) 
    !>  (default: `0` = smooth function)
    integer, optional, intent(in) :: derivative
    !> weighting function   
    !> valid options: `'rectangle'` (no weighting)   
    !>                `'hann'` (Hann function, default)
    character(len=*), optional, intent(in) :: weighting 
  
    character(len=*), parameter :: PROCEDURE_NAME = 'savgol_single'

    integer :: dord
    real(kind=dp), pointer :: y_many(:,:)

    dord = 0
    if (present(derivative)) dord = derivative  
    y_many(1:size(y), 1:1) => y

    call savgol_many( x, y_many, degree=degree, window=window, derivative=[dord], weighting=weighting )
  end subroutine savgol_single

  !> Savitzky-Golay filter for equidistantly sampled function. See [[savgol_many(subroutine)]].
  subroutine savgol_equi( x, y, deg, win, dord, wfun )
    !> values \(x_i\) in increasing order
    real(kind=dp), intent(in) :: x(:)
    !> on input: values \(y_i\) in first column   
    !> on output: smooth approximation \(f(x_i)\) or selected derivatives \(f^{(d)}(x_i)\)
    real(kind=dp), intent(inout) :: y(:,:)
    !> polynomial degree 
    integer, intent(in) :: deg
    !> window length \(m\) (\(2m+1\) points centered around \(x_i\) are considered for least squares fit) 
    integer, intent(in) :: win
    !> order of derivatives to compute (derivatives of order greater than polynomial degree are zero) 
    integer, intent(in) :: dord(:)
    !> weighting function   
    character(len=*), intent(in) :: wfun
  
    character(len=*), parameter :: PROCEDURE_NAME = 'savgol_equi'

    integer :: n, d, i
    real(kind=dp) :: linfun(0:1)
    real(kind=dp) :: dx

    real(kind=dp), allocatable :: y0(:)
    real(kind=dp), allocatable, target :: w(:,:)
    real(kind=dp), pointer :: wgt(:,:)

    n = size(x)

    ! make copy of input data
    y0 = y(1:n, 1)

    ! subtract linear function to turn endpoints to zero
    linfun(0) = (x(n)*y0(1) - x(1)*y0(n)) / (x(n) - x(1))
    linfun(1) = (y0(n) - y0(1)) / (x(n) - x(1))
    y0 = y0 - (linfun(0) + linfun(1)*x)

    ! get weights
    w = equi_weights( deg, win, wfun )
    wgt(-win:win, 0:deg) => w

    ! compute smooth function and derivatives
    do d = 1, size(dord)
      dx = (x(n) - x(1)) / (n - 1)
      ! inner points
      do i = win + 1, n - win
        y(i, d) = dot_product( wgt(:, dord(d)), y0(i-win:i+win) )
      end do
      ! boundary points
      do i = win, 1, -1
        y(i, d) = dot_product( wgt(:, dord(d)), [-y0(2+win-i:2:-1), y0(1:win+i)] )
      end do
      do i = n - win + 1, n
        y(i, d) = dot_product( wgt(:, dord(d)), [y0(i-win:n), -y0(n-1:2*n-i-win:-1)] )
      end do
      ! scale for derivatives
      if (dord(d) > 0) y(1:n, d) = y(1:n, d) * (1.0_dp / dx**dord(d)) 
      ! add back linear function
      if (dord(d) == 0) then
        y(1:n, d) = y(1:n, d) + (linfun(0) + linfun(1)*x) 
      else if (dord(d) == 1) then
        y(1:n, d) = y(1:n, d) + linfun(1) 
      end if
    end do
  end subroutine savgol_equi

  !> Interpolation weigths \(w_{jk}\) for equidistant sampling, given polynomial degree \(p\) and 
  !> window length \(m\) such that the coeficients \(a_{ki}\) of the fitted polynomial for 
  !> the points centered around \(x_i\) are given by
  !> \[ a_{ki} = \sum_{j=-m}^{m} w_{jk}\, y_{i+j} \;. \]
  !> The localy fitted polynomial is then given by
  !> \[ q(x) = \sum_{k=0}^p a_{ki} \, \left(\frac{x-x_i}{\Delta x}\right)^k \;, \]
  !> where \(\Delta x\) is the spacing between the points \(x_i\).
  function equi_weights( degree, window, weighting ) result( weights )
    !> polynomial degree 
    integer, intent(in) :: degree
    !> window length \(m\) (\(2m+1\) points centered around \(x_i\) are considered for least squares fit) 
    integer, intent(in) :: window
    !> weighting function   
    !> valid options: `'rectangle'` (no weighting)   
    !>                `'hann'` (Hann function)
    character(len=*), intent(in) :: weighting 
    !> interpolation weights \(w_{jk}\)
    real(kind=dp), allocatable :: weights(:,:)
  
    character(len=*), parameter :: PROCEDURE_NAME = 'equi_weights'

    integer :: j, k, info

    integer, allocatable :: ipiv(:)
    real(kind=dp), allocatable :: lhs(:,:), rhs(:,:)

    ! right hand side matrix
    allocate( rhs(0:degree, -window:window), lhs(0:degree, 0:degree) )
    rhs(0, :) = 1.0_dp
    if (degree > 0) rhs(1, :) = [(real(j, kind=dp), j=-window, window)]
    do k = 2, degree
      rhs(k, :) = rhs(k-1, :) * rhs(1, :)
    end do
    call apply_weighting_function
    ! left hand side matrix
    call dgemm( 'n', 't', degree+1, degree+1, 2*window+1, 1.0_dp, &
      rhs, degree+1, &
      rhs, degree+1, 0.0_dp, &
      lhs, degree+1 )
    call apply_weighting_function
    ! solve linear system for weights
    allocate( ipiv(2*window+1) )
    call dgesv( degree+1, 2*window+1, lhs, degree+1, ipiv, rhs, degree+1, info )
    if (info /= 0) error stop 'Error ('//MODULE_NAME//':'//PROCEDURE_NAME//'): '// &
      'Failed to solve linear system for interpolation weights.'
    allocate( weights(-window:window, 0:degree) )
    weights = transpose( rhs )

  contains

    subroutine apply_weighting_function
      use constants, only: pi
      integer :: j
      if (weighting == 'hann') then
        do j = -window, window
          rhs(:, j) = rhs(:, j) * cos( 0.5_dp * pi * real(j, kind=dp) / real(window, kind=dp) )
        end do
      end if
    end subroutine apply_weighting_function
  end function equi_weights
end module savitzky_golay
