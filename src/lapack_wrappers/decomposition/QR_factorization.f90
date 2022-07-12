!> Interfaces for QR factorization. The interfaces combine wrappers for the LAPACK routines
!> DGEQP3, ZGEQP3
module qr_factorization
  use asserts, only: assert
  use modmpi, only: terminate_if_false
  use precision, only: dp
  use constants, only: zzero
  use lapack_f95_interfaces, only: dgeqp3, zgeqp3, dorgqr, zungqr

  implicit none
  
  private
  public :: qr_column_pivot, xgeqp3, xorgqr, extract_R


  !> Calculate the QR factorization with column pivoting 
  !> of a complex matrix \( A \in \mathcal{K}^{m \times n} \) such that
  !> \[
  !>   A \cdot P = Q \cdot R,
  !> \]
  !> where P is a column permutation matrix, \( Q \in \mathcal{K}^{m \times m} \) is a 
  !> unitary matrix and \( R \in \mathcal{K}^{l \times k} \) is a upper triangular 
  !> matrix, where \( k = \min(m, n) \) and \( l = \max(m, n) \).
  !> If 'm < n', the QR factorization is not defined uniquely and therefore not supported.
  !> \( \mathcal{K} \) is one of \( \mathcal{R} \) or \( \mathcal{C} \).
  interface qr_column_pivot
    module procedure qr_column_pivot_real_dp, &
                     qr_column_pivot_complex_dp
  end interface qr_column_pivot

  !> Calculate the QR factorization with column pivoting. See [[qr_column_pivot_real_dp]].
  !> Use this routine if \( \mathbf{A} \) can be mutated. It will be overwritten with \( \mathbf{R} \)
  !> and the reflector representation of \( \mathbf{Q} \). To get the actual matrix, refer to
  !> [[dorgqr_wrapper]].
  !>
  !> This routine acts on the arrays as expected by [[*geqp3]].
  interface xgeqp3
    module procedure dgeqp3_wrapper, &
                     zgeqp3_wrapper
  end interface xgeqp3
  
  !> Calculate the orthonormal matrix \( \mathbf{Q} \) from the reflector representation
  !> (as reurned by [[dgeqp3_wrapper]] and
  !> [[zgeqp3_wrapper]] respectively), which is defined as the
  !> first 'n' columns of a product of 'k' elementary reflectors of order 'm':
  !> \[
  !>    \mathbf{Q} = \mathbf{H}(1) \mathbf{H}(2) ...  \mathbf{H}(k).
  !> \]
  !>
  !> This routine acts on the arrays as expected by [[*orgqr]].
  interface xorgqr
    module procedure dorgqr_wrapper, &
                     zorgqr_wrapper
  end interface xorgqr


  !> Extract from a matrix \( \mathbf{A} \in \mathcal{K}^{m \times n} \) the 
  !> upper trapezoidal matrix \( \mathbf{R} \in \mathcal{R}^{\min(m, n), n})
  !> two arrays \( \mathbf{L} \) \( \mathbf{U} \).
  interface extract_R
    module procedure extract_R_real_dp, &
                     extract_R_complex_dp
  end interface extract_R

contains

  !> Calculate the QR factorization with column pivoting 
  !> of a real matrix \( \mathbf{A} \in \mathcal{R}^{m \times n} \), where 'm \ge n', such that
  !> \[
  !>   \mathbf{A} \cdot \mathbf{P} = \mathbf{Q} \cdot \mathbf{R},
  !> \]
  !> where \( \mathbf{P} \) is a column permutation matrix, \( \mathbf{Q} \in \mathcal{R}^{m \times m} \) is a 
  !> unitary matrix and \( \mathbf{R} \in \mathcal{R}^{l \times k} \) is a upper triangular 
  !> matrix, where \( k = \min(m, n) \) and \( l = \max(m, n) \).
  !> If 'm < n', the QR factorization is not defined uniquely and therefore not supported.
  !> Use this interface if \( \mathbf{A} \) must not be mutated.
  subroutine qr_column_pivot_real_dp(A, P, Q, R)
    !> Input matrix \(\mathbf{A}\)
    real(dp), contiguous, intent(in) :: A(:, :)
    !> Permutation map corresponding to permutation matrix \(\mathbf{P}\)
    integer, contiguous, intent(out) :: P(:)
    !> Unitary matrix \( \mathbf{Q} \), this will only be calculated if 'm > n'.
    real(dp), intent(out), contiguous, optional :: Q(:, :) 
    !> Trapezoidal matrix \( \mathbf{R} \)
    real(dp), intent(out), contiguous, optional :: R(:, :)

    integer :: m, n
    real(dp), allocatable :: A_(:, :), tau(:)

    m = size(A, dim=1); n = size(A, dim=2)

! This preprocessor usage is deliberate, to prevent if statements being evaluated in production code.
#ifdef USE_ASSERT
    call assert(m >= n, 'A has m < n.')
    if (present(R)) then
      call assert(all(shape(R) == [m, min(m, n)]), 'R must be of shape [m, min(m, n)].')
    end if 

    if (present(Q)) then
      call assert(all(shape(Q) == [m, m]), 'Q must be a square matrix of size m.')
    end if
#endif

    A_ = A; allocate(tau(min(m, n)))
    P = 0
    call xgeqp3(A_, P, tau)

    if (present(R)) call extract_R(A_, R)
    if (present(Q)) then 
      Q(:, 1 : n) = A_
      if(m  >n) Q(:, n + 1 : m) = 0.0_dp
      call xorgqr(Q, tau)
    end if
  end subroutine qr_column_pivot_real_dp

  !> Calculate the QR factorization with column pivoting. See [[qr_column_pivot_real_dp]].
  !> Use this routine if \( \mathbf{A} \) can be mutated. It will be overwritten with \( \mathbf{R} \)
  !> and the reflector representation of \( \mathbf{Q} \). To get the actual matrix, refer to 
  !> [[dorgqr_wrapper]].
  !>
  !> This routine acts on the arrays as expected by [[dgeqp3]].
  subroutine dgeqp3_wrapper(A, jpvt, tau)
    !> On entry input matrix \( \mathbf{A} \), on exit the upper triangle contains 
    !> the \(\min(m, n) \times n\) upper trapexzoidal matrix \( \mathbf{R} \).
    !> The elements below the diagonal, together with **tau**, represent the unitary 
    !> matrix \( \mathbf{Q} \) as a product of \( \min(m, n) \) reflectors.
    real(dp), contiguous, intent(inout) :: A(:, :)
    !> On entry, if **jpvt**\((j) \ne 0 \), the 'j'-th column of \( \mathbf{A}\)
    !> is permuted to the front. If **jpvt**\((j) = 0 \), the 'j'-th column of 
    !> \( \mathbf{A}\) is a free column.
    !> On exit, **jpvt** contains the permutation map corresponding to permution 
    !> matrix \(mathbf{P}\).
    integer, contiguous, intent(inout) :: jpvt(:)
    !> Scalar factors of the elementary reflectors
    real(dp), contiguous, intent(out) :: tau(:)

    integer :: m, n, lwork, info
    real(dp), allocatable :: work(:)

    m = size(A, dim=1); n = size(A, dim=2)

    call assert(m >= n, 'A has m < n.')
    call assert(size(jpvt) == n, &
               'The number of elements of jpvt is not equal to the number of culumns of A.')
    call assert(size(tau) == min(m, n), &
               'The number of elements of taus is not equal to the mininum of number of rows and culomns of A.')

    lwork = -1; allocate(work(1))
    call dgeqp3(m, n, A, m, jpvt, tau, work, lwork, info)
    call terminate_if_false(info == 0, 'dgqp3 work space query failed.')

    lwork = work(1); deallocate(work); allocate(work(lwork))
    call dgeqp3(m, n, A, m, jpvt, tau, work, lwork, info)
    call terminate_if_false(info == 0, 'dgqp3 QR factorization failed.')
  end subroutine dgeqp3_wrapper

  !> Calculate the orthonormal matrix \( \mathbf{Q} \) from the reflector representation
  !> as reurned by [[dgeqp3_wrapper]], which is defined as the
  !> first 'n' columns of a product of 'k' elementary reflectors of order 'm':
  !> \[
  !>    \mathbf{Q} = \mathbf{H}(1) \mathbf{H}(2) ...  \mathbf{H}(k).
  !> \]
  subroutine dorgqr_wrapper(A, tau)
    !> On entry the 'i'-th column must contain the the elementary reflector \( \mathbf{H}(i) \),
    !> as returned by [[dgeqp3_wrapper]]. On exit, the array
    !> contains \( \mathbf{Q} \)
    real(dp), contiguous, intent(inout) :: A(:, :)
    !> Scalar factors of the elementary reflectors
    real(dp), contiguous, intent(in) :: tau(:)

    integer :: m, n, k, lwork, info
    real(dp), allocatable :: work(:)

    m = size(A, dim=1); n = size(A, dim=2); k = size(tau)

    call assert(m >= n, 'A has m < n.')

    lwork = -1; allocate(work(1))
    call dorgqr(m, n, k, A, m, tau, work, lwork, info)
    call terminate_if_false(info == 0, 'dorgqr work space query failed.')
    
    lwork = work(1); deallocate(work); allocate(work(lwork))
    call dorgqr(m, n, k, A, m, tau, work, lwork, info)
    call terminate_if_false(info == 0, 'dorgqr calculating Q from elementary reflector failed.')
  end subroutine dorgqr_wrapper


  !> Calculate the QR factorization with column pivoting 
  !> of a complex matrix \( \mathbf{A} \in \mathcal{C}^{m \times n} \) such that
  !> \[
  !>   \mathbf{A} \cdot \mathbf{P} = \mathbf{Q} \cdot \mathbf{R},
  !> \]
  !> where \( \mathbf{P} \) is a column permutation matrix, \( \mathbf{Q} \in \mathcal{C}^{m \times m} \) is a 
  !> unitary matrix and \( \mathbf{R} \in \mathcal{C}^{l \times k} \) is a upper triangular 
  !> matrix, where \( k = \min(m, n) \) and \( l = \max(m, n) \).
  !> If 'm < n', the QR factorization is not defined uniquely and therefore not supported.
  !> Use this interface if \( \mathbf{A} \) must not be mutated.
  subroutine qr_column_pivot_complex_dp(A, P, Q, R)
    !> Input matrix \(\mathbf{A}\)
    complex(dp), contiguous, intent(in) :: A(:, :)
    !> Permutation map corresponding to permutation matrix \(\mathbf{P}\)
    integer, contiguous, intent(out) :: P(:)
    !> Unitary matrix \( \mathbf{Q} \)
    complex(dp), intent(out), contiguous, optional :: Q(:, :) 
    !> Trapezoidal matrix \( \mathbf{R} \)
    complex(dp), intent(out), contiguous, optional :: R(:, :)

    integer :: m, n
    complex(dp), allocatable :: A_(:, :), tau(:)

    m = size(A, dim=1); n = size(A, dim=2)

! This preprocessor usage is deliberate, to prevent if statements being evaluated in production code.
#ifdef USE_ASSERT
    if (present(R)) then
      call assert(all(shape(R) == [m, min(m, n)]), 'R must be of shape [m, min(m, n)].')
    end if 

    if (present(Q) .and. m > n) then
      call assert(all(shape(Q) == [m, m]), 'Q must be a square matrix of size m.')
    end if
#endif

    A_ = A; allocate(tau(min(m, n)))
    P = 0
    call xgeqp3(A_, P, tau)

    if (present(R)) call extract_R(A_, R)
    if (present(Q)) then 
      Q(:, 1 : n) = A_
      if(m>n) Q(:, n+1 : m) = zzero
      call xorgqr(Q, tau)
    end if 
  end subroutine qr_column_pivot_complex_dp

  !> Calculate the QR factorization with column pivoting. See [[qr_column_pivot_real_dp]].
  !> Use this routine if \( \mathbf{A} \) can be mutated. It will be overwritten with \( \mathbf{R} \)
  !> and the reflector representation of \( \mathbf{Q} \). To get the actual matrix, refer to 
  !> [[zorgqr_wrapper]].
  !>
  !> This routine acts on the arrays as expected by [[zgeqp3]].
  subroutine zgeqp3_wrapper(A, jpvt, tau)
    !> On entry input matrix \( \mathbf{A} \), on exit the upper triangle contains 
    !> the \(\min(m, n) \times n\) upper trapexzoidal matrix \( \mathbf{R} \).
    !> The elements below the diagonal, together with **tau**, represent the unitary 
    !> matrix \( \mathbf{Q} \) as a product of \( \min(m, n) \) reflectors.
    complex(dp), contiguous, intent(inout) :: A(:, :)
    !> On entry, if **jpvt**\((j) \ne 0 \), the 'j'-th column of \( \mathbf{A}\)
    !> is permuted to the front. If **jpvt**\((j) = 0 \), the 'j'-th column of 
    !> \( \mathbf{A}\) is a free column.
    !> On exit, **jpvt** contains the permutation map corresponding to permution 
    !> matrix \(mathbf{P}\).
    integer, contiguous, intent(inout) :: jpvt(:)
    !> Scalar factors of the elementary reflectors
    complex(dp), contiguous, intent(out) :: tau(:)

    integer :: m, n, lwork, info
    real(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: work(:)

    m = size(A, dim=1); n = size(A, dim=2)

    !call assert(m >= n, 'A has m < n.')
    call assert(size(jpvt) == n, &
               'The number of elements of jpvt is not equal to the number of culumns of A.')
    call assert(size(tau) == min(m, n), &
               'The number of elements of taus is not equal to the mininum of number of rows and culomns of A.')

    lwork = -1; allocate(work(1)); allocate(rwork(2 * n))
    call zgeqp3(m, n, A, m, jpvt, tau, work, lwork, rwork, info)
    call terminate_if_false(info == 0, 'zgqp3 work space query failed.')

    lwork = work(1); deallocate(work); allocate(work(lwork))
    call zgeqp3(m, n, A, m, jpvt, tau, work, lwork, rwork, info)
    call terminate_if_false(info == 0, 'zgqp3 QR factorization failed.')
  end subroutine zgeqp3_wrapper

  !> Calculate the orthonormal matrix \( \mathbf{Q} \) from the reflector representation
  !> as reurned by [[dgeqp3_wrapper]], which is defined as the
  !> first 'n' columns of a product of 'k' elementary reflectors of order 'm':
  !> \[
  !>    \mathbf{Q} = \mathbf{H}(1) \mathbf{H}(2) \ddots  \mathbf{H}(k).
  !> \]
  subroutine zorgqr_wrapper(A, tau)
    !> On entry the 'i'-th column must contain the the elementary reflector \( \mathbf{H}(i) \),
    !> as returned by [[dgeqp3_wrapper]]. On exit, the array
    !> contains \( \mathbf{Q} \)
    complex(dp), contiguous, intent(inout) :: A(:, :)
    !> Scalar factors of the elementary reflectors
    complex(dp), contiguous, intent(in) :: tau(:)

    integer :: m, n, k, lwork, info
    complex(dp), allocatable :: work(:)

    m = size(A, dim=1); n = size(A, dim=2); k = size(tau)

    call assert(m >= n, 'A has m < n.')

    lwork = -1; allocate(work(1))
    call zungqr(m, n, k, A, m, tau, work, lwork, info)
    call terminate_if_false(info == 0, 'dorgqr work space query failed.')
    
    lwork = work(1); deallocate(work); allocate(work(lwork))
    call zungqr(m, n, k, A, m, tau, work, lwork, info)
    call terminate_if_false(info == 0, 'dorgqr calculating Q from elementary reflector failed.')
  end subroutine zorgqr_wrapper


! utils

  !> Extract from a matrix \( \mathbf{A} \in \mathcal{R}^{m \times n} \) the 
  !> upper trapezoidal matrix \( \mathbf{R} \in \mathcal{R}^{\min(m, n), n}).
  subroutine extract_R_real_dp(A, R) 
    !> Input matrix \( \mathbf{A} \); contains 
    !> the factors \( \mathbf{L} \) and \( \mathbf{U} \), where the 
    !> unit diagonals of \( \mathbf{L} \) are not stored.
    real(dp), intent(in), contiguous :: A(:, :)
    !> Lower triangular matrix \( \mathbf{L} \)
    real(dp), intent(out), contiguous :: R(:, :)

    integer :: i, m, n, k, l

    m = size(A, dim=1)
    n = size(A, dim=2)
    k = min(m, n)
    l = max(m, n)

    call assert(m >= n, 'A has m < n.')
    call assert(all(shape(R) == [max(m, n), min(m, n)]), &
               'R must be of shape [max(m, n), min(m, n)].')

    do i=1, k-1
      R(1 : i, i) = A(1 : i, i)
      R(i+1 : l, i) = 0.0_dp
    end do

    R(1 : k, k) = A(1 : k, k)

    if (m > n) then
      R(k+1 : m, :) = 0.0_dp
    end if 
  end subroutine extract_R_real_dp

  !> Extract from a matrix \( \mathbf{A} \in \mathcal{C}^{m \times n} \) the 
  !> upper trapezoidal matrix \( \mathbf{R} \in \mathcal{C}^{\min(m, n), n}).
  subroutine extract_R_complex_dp(A, R) 
    !> Input matrix \( \mathbf{A} \); contains 
    !> the factors \( \mathbf{L} \) and \( \mathbf{U} \), where the 
    !> unit diagonals of \( \mathbf{L} \) are not stored.
    complex(dp), intent(in), contiguous :: A(:, :)
    !> Lower triangular matrix \( \mathbf{L} \)
    complex(dp), intent(out), contiguous :: R(:, :)

    integer :: i, m, n, k, l

    m = size(A, dim=1)
    n = size(A, dim=2)
    k = min(m, n)
    l = max(m, n)

    call assert(all(shape(R) == [max(m, n), min(m, n)]), &
               'R must be of shape [max(m, n), min(m, n)].')

    do i=1, k-1
      R(1 : i, i) = A(1 : i, i)
      R(i+1 : l, i) = zzero
    end do

    R(1 : k, k) = A(1 : k, k)

    if (m > n) then
      R(k+1 : m, :) = zzero
    end if
  end subroutine extract_R_complex_dp

end module qr_factorization