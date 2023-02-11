!> Interfaces for singular value decomposition (svd) and matrix rank calculation, 
!> relying on svd. The interface for svd combines wrappers for the LAPACK routines
!> **[[dgesdd]]**, **[[zgesdd]]**.
module singular_value_decomposition
  use precision, only: dp
  use asserts, only: assert
  use modmpi, only: mpiglobal, terminate_if_false
  use lapack_f95_interfaces, only: dgesdd, zgesdd

  implicit none 
  
  private
  public :: svd_divide_conquer, xgesdd

  !> Default tolerance \( \frac {\min(m, n)} {\max(m, n)} < \text{tol} \) such that 
  !> \( {\min(m, n)}  \ll {\max(m, n)} \) can be assumed
  real(dp), parameter :: default_tol_narrow_matrix = 1e-5_dp

  !> Calculate the singular value decomposition (SVD) of a matrix \( \mathbf{A} \):
  !> \[
  !>    \mathbf{A} = \mathbf{U} \cdot \Sigma \cdot \mathbf{V}^\dagger
  !> \]
  !> where \( \Sigma \) is an m-by-n matrix which is zero except for its
  !> \( \min(m,n) \) diagonal elements, \( \mathbf{U} \) is an \( m \)-by-\( m \) unitary matrix, and
  !> \( \mathbf{V} \) is an \( n \)-by-\( n \) unitary matrix.  The diagonal elements of \( \Sigma \)
  !> are the singular values of \( \mathbf{A} \); they are real and non-negative, and
  !> are returned in descending order as an array.  The first \( \min(m,n) \) columns of
  !> \( \mathbf{U} \) and \( \mathbf{V} \) are the left and right singular vectors of \( \mathbf{A} \)
  !> respectively.
  !>
  !> Note that the routine returns \( \mathbf{V}^\top \), not \( \mathbf{V} \).
  interface svd_divide_conquer
    module procedure svd_divide_and_conquer_real_dp, &
                     svd_divide_and_conquer_complex_dp
  end interface svd_divide_conquer


  !> See [[svd_divide_and_conquer]]. Use this interface, if \(\mathbf{A}\) can be mutated.
  !>
  !> This routine acts on the arrays as expected by [[*gesdd]].
  interface xgesdd
    module procedure dgesdd_wrapper, &
                     zgesdd_wrapper
  end interface xgesdd

  contains 

! singular value decomposition with divide and conquer algorithm

  !> Calculate the singular value decomposition (SVD) of a matrix \( \mathbf{A} \in \mathbb{R}^{m \times n} \).
  !> The SVD is written as
  !> \[
  !>    \mathbf{A} = \mathbf{U} \cdot \Sigma \cdot \mathbf{V}^\top
  !> \]
  !> where \( \Sigma \in \mathbb{R}^{m \times n} \) with \(\sigma_{ij} = 0 \) for \( i \neq j \), 
  !> \( \mathbf{U} \in \mathbb{R}^{m \times m} \) and \( \mathbf{V} \in \mathbb{R}^{n \times n} \) are  orthogonal matrices.  
  !> The diagonal elements of \( \Sigma \) are the singular values of \( \mathbf{A} \); they are real and non-negative, and
  !> are returned in descending order.  The first \( \min(m,n) \) columns of
  !> \( \mathbf{U} \) and \( \mathbf{V} \) are the left and right singular vectors of \( \mathbf{A} \).
  !>
  !> Note that the routine returns \( \mathbf{V}^\top \), not \( \mathbf{V} \).
  !>
  !> The divide and conquer algorithm makes very mild assumptions about
  !> floating point arithmetic. It will work on machines with a guard
  !> digit in add/subtract, or on those binary machines without guard
  !> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
  !> Cray-2. It could conceivably fail on hexadecimal or decimal machines
  !> without guard digits, but we know of none.
  !> Use this interface, if \(\mathbf{A}\) should be unchanged.
  !> (This documation is from netlib.org)
  subroutine svd_divide_and_conquer_real_dp(A, sigma, U, V_T)
    !> Input matrix \( \mathbf{A} \in \mathbb{R}^{m \times n} \)
    real(dp), intent(in), contiguous :: A(:, :)
    !> Singular values (diagonal elements of \( \Sigma \)), stored as array, 
    !> sorted such that \(\sigma_i \ge \sigma_{i+1}\)
    real(dp), intent(out), contiguous :: sigma(:)
    !> On exit, contains (part of the) left singular vectors \( \mathbf{U} \in \mathbb{R}^{m \times m} \).
    !> The vectors are stored column-wise in the array **U**. 
    !>
    !> The shape of **U** determines how many left singular vectors are calculated.
    !>
    !>  - If **U** is not given, only **sigma** is calculated:
    !> 
    !>  - If \(\text{shape} ( \)**U** \( ) =   [m , m] \), 
    !>    all left singular vectors are calculated
    !>
    !>  - If \( \text{shape} ( \)**U** \( ) =   [m , \min(m, n)] \),
    !>    the first \( \min(m, n) \) left singular vectors are calculated.
    real(dp), intent(out), contiguous, optional :: U(:,:)
    !> On exit, contains (part of the) transpose of the right singular vectors \( \mathbf{V}^\top \).
    !> The vectors are stored row-wise in the array **V_T**.
    !>
    !> The shape of **V_T** determines how many the right singular vectors are calculated.
    !>
    !> - If **V_T** is not given, only **sigma** is calculated:
    !>
    !> - If \(\text{shape} ( \)**V_T** \( ) =   [n , n] \),
    !>   all right singular vectors are calculated
    !>
    !> - If \(m < n \) and \(\text{shape} ( \)**V_T** \( ) =   [\min(m, n) , n] \),
    !>   the first \( \min(m, n) \) right singular vectors are calculated.
    real(dp), intent(out), contiguous, optional :: V_T(:,:)
        
    character(1) :: jobz
    real(dp), allocatable :: A_(:, :), U_(:, :), V_T_(:, :)

    A_ = A

    if(present(U) .and. present(V_T)) then
      jobz = setup_jobz(shape(A), shape(U), shape(V_T))
      call dgesdd_wrapper(jobz, A_, sigma, U, V_T)
    
    else if(.not. (present(U) .or. present(V_T))) then
      jobz = 'N'
      allocate(U_(1, 1)); allocate(V_T_(1, 1))
      call dgesdd_wrapper(jobz, A_, sigma, U_, V_T_)
    
    else 
      call terminate_if_false(.false., 'Either both, U and V_T, must be present or none of both.')
    end if

    if (jobz == 'O' .and. size(A, dim=1) >= size(A, dim=2)) U = A_ 
    if (jobz == 'O' .and. size(A, dim=1)  < size(A, dim=2)) V_T = A_ 
  end subroutine svd_divide_and_conquer_real_dp


  !> See [[svd_divide_and_conquer_real_dp]]. Use this interface if \(\mathbf{A}\) can be mutated.
  !>
  !> This routine acts on the arrays as expected by [[dgesdd]].
  subroutine dgesdd_wrapper(jobz, A, sigma, U, V_T)
    !> Specifies options for computing all or part of the matrix U (see zgesdd):
    !>
    !> - `'A'`: All \( m \) left and all \( n \) right singular values are calculated.
    !>   The left singular vectors are stored column-wise in **U**.
    !>   The complex-conjugate of the right singular vectors are stored row-wise in **V_T**.
    !>
    !> - `'S'`: The first \( \min(m, n) \) clumns of **U** and rows of **V_T** are calculated
    !>
    !> - `'O'`:
    !>
    !>     - If \( m \ge n \):
    !>
    !>         - The first \( n \) columns of **A** are overwritten with  the first \( n \)
    !>           left singular vectors.
    !>
    !>         - **U** is not referenced.
    !>
    !>         - **V_T** contains all right singular vectors, stored row-wise.
    !>
    !>     - If \( m < n \):
    !>
    !>         - The first \( m \) rows of **A** are overwritten with the
    !>           the first \( m \) right singular vectors.
    !>
    !>         - **U** contains all left singular vectors, stored column-wise.
    !>
    !>         - **V_T** is not referenced.
    !>
    !> - `'N'`:  Neither left nor right singular vectors are calculated.
    character :: jobz
    !> On entry, contains the matrix \( \mathbf{A} \in \mathbb{R}^{m \times n}\).
    !>
    !> On exit:
    !>
    !> - If `jobz == 'O'`:
    !>
    !>     - If \( m \ge n \), the first \( n \) columns of **A** are overwritten with  the first \( n \)
    !>       left singular vectors.
    !>
    !>     - If \( m < n \), the first \( m \) rows of **A** are overwritten with the
    !>       the first \( m \) right singular vectors.
    !>
    !> - If `jobz \= 'O'`, the contents of A are destroyed.
    real(dp), intent(inout), contiguous :: A(:, :)
    !> Singular values (diagonal elements of \( \Sigma \)), stored as array, sorted such that \(\sigma_i \ge \sigma_{i+1}\)
    real(dp), intent(out), contiguous :: sigma(:)
    !> Contains (part of) the left singular vectors \( \mathbf{U} \in \mathbb{R}^{m \times m} \).
    !> The vectors are stored column-wise in the array **U**.
    !>
    !> - If `jobz == 'A'`, **U** contains all left singular vectors, stored column-wise.
    !>
    !> - If `jobz == 'S'`, **U** contains the first \( \min(m, n) \) left singular vectors, stored column-wise.
    !>
    !> - If `jobz == 'O'`:
    !>
    !>     - If \( m \ge n \), **U** is not referenced.
    !>
    !>     - If \( m < n \), **U** contains all left singular vectors, stored column-wise.
    !>
    !> - If `jobz == 'N'`, **U** is not referenced.
    real(dp), intent(out), contiguous :: U(:,:)
    !> Contains (part of) the transpose of the right singular vectors \( \mathbf{V}^\top \in \mathbb{R}^{n \times n} \).
    !> The vectors are stored row-wise in the array **V_T**.
    !>
    !> - If `jobz == 'A'`, **V_T** contains all left singular vectors, stored column-wise.
    !>
    !> - If `jobz == 'S'`, **V_T** contains the first \( \min(m, n) \) left singular vectors, stored column-wise.
    !>
    !> - If `jobz == 'O'`:
    !>
    !>     - If \( m \ge n \), **V_T** contains the all right singular vectors, stored row-wise.
    !>
    !>     - If \( m < n \), **V_T** is not referenced.
    !>
    !> - If `jobz == 'N'`, **V_T** is not referenced.
    real(dp), intent(out), contiguous :: V_T(:,:)
    
    integer :: m, n, k, info, lwork
    integer, allocatable :: iwork(:)
    real(dp), allocatable :: work(:)


    m = size(A, dim=1)
    n = size(A, dim=2)
    k = min(m, n)

    call assert_array_shape(jobz, m, n, k, size(sigma), shape(U), shape(V_T))
    
    lwork = -1
    allocate(work(1))
    allocate(iwork(8 * k))
      
    call dgesdd(jobz, m, n, A, m, sigma, U, size(U, dim=1), V_T, size(V_T, dim=1), work, lwork, iwork, info)
    call terminate_if_false(info == 0, &
                           'dgesdd failed for searching optiomal work space.')

    lwork = idnint(work(1))
    deallocate(work); allocate(work(lwork))
    
    call dgesdd(jobz, m, n, A, m, sigma, U, size(U, dim=1), V_T, size(V_T, dim=1), work, lwork, iwork, info)

    call terminate_if_false(info == 0, &
                           'dgesdd failed for calculating the SVD of A.')
  end subroutine dgesdd_wrapper


  !> Calculate the singular value decomposition (SVD) of a matrix 
  !> \( \mathbf{A} \in \mathbb{C}^{m \times n} \).
  !> The SVD is written as
  !> \[
  !>    \mathbf{A} = \mathbf{U} \cdot \Sigma \cdot V^\dagger
  !> \]
  !> where \( \Sigma \in \mathbb{R}^{m \times n} \) with \(\sigma_{ij} = 0 \) for
  !> \( i \neq j \), \( \mathbf{U} \in \mathbb{C}^{m \times m} \) and
  !> \( \mathbf{V} \in \mathbb{C}^{n \times n} \) are  orthogonal matrices.  The diagonal elements of
  !> \( \Sigma \) are the singular values of \( \mathbf{A} \); they are real and non-negative, and
  !> are returned in descending order.  The first \( \min(m,n) \) columns of
  !> \( \mathbf{U} \) and \( \mathbf{V} \) are the left and right singular vectors of \( \mathbf{A} \).
  !>
  !> Note that the routine returns \( \mathbf{V}^\dagger \), not \( \mathbf{V} \).
  !>
  !> The divide and conquer algorithm makes very mild assumptions about
  !> floating point arithmetic. It will work on machines with a guard
  !> digit in add/subtract, or on those binary machines without guard
  !> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
  !> Cray-2. It could conceivably fail on hexadecimal or decimal machines
  !> without guard digits, but we know of none.
  !> Use this interface, if \(\mathbf{A}\) should be unchanged.
  !> (This documation is from netlib.org)
  subroutine svd_divide_and_conquer_complex_dp(A, sigma, U, V_H)
    !> Input matrix \( \mathbf{A} \in \mathbb{C}^{m \times n} \)
    complex(dp), intent(in), contiguous :: A(:,:)
    !> Singular values (diagonal elements of \( \Sigma \)), stored as array,
    !> sorted such that \(\sigma_i \ge \sigma_{i+1}\)
    real(dp), intent(out), contiguous :: sigma(:)
    !> Contains (part of) the left singular vectors \( \mathbf{U} \in \mathbb{C}^{m \times m} \).
    !> The vectors are stored column-wise in the array **U**.
    !>
    !> The shape of **U** determines how many left singular vectors are calculated.
    !>
    !>  - If **U** is not given, only **sigma** is calculated:
    !>
    !>  - If \(\text{shape} ( \)**U** \( ) =   [m , m] \),
    !>    all left singular vectors are calculated
    !>
    !>  - If \( \text{shape} ( \)**U** \( ) =   [m , \min(m, n)] \),
    !>    the first \( \min(m, n) \) left singular vectors are calculated.
    complex(dp), intent(out), contiguous, optional :: U(:,:)
    !> Contains (part of) the transpose-conjugate of the right singular vectors \( \mathbf{V}^\dagger \in \mathbb{C}^{n \times n} \).
    !> The conjugated vectors are stored row-wise in the array **V_H**.
    !>
    !> The shape of **V_H** determines how many the right singular vectors are calculated.
    !>
    !> - If **V_H** is not given, only **sigma** is calculated:
    !>
    !> - If \(\text{shape} ( \)**V_H** \( ) =   [n , n] \),
    !>   all right singular vectors are calculated
    !>
    !> - If \(m < n \) and \(\text{shape} ( \)**V_H** \( ) =   [\min(m, n) , n] \),
    !>   the first \( \min(m, n) \) right singular vectors are calculated.
    complex(dp), intent(out), contiguous, optional :: V_H(:,:)
        
    character(1) :: jobz
    complex(dp), allocatable :: A_(:, :), U_(:, :), V_H_(:, :)

    A_ = A

    if(present(U) .and. present(V_H)) then
      jobz = setup_jobz(shape(A), shape(U), shape(V_H))
      call zgesdd_wrapper(jobz, A_, sigma, U, V_H)
    else
      jobz = 'N'
      allocate(U_(1, 1)); allocate(V_H_(1, 1))
      call zgesdd_wrapper(jobz, A_, sigma, U_, V_H_)
    end if

    if (jobz == 'O' .and. size(A, dim=1) >= size(A, dim=2)) U = A_ 
    if (jobz == 'O' .and. size(A, dim=1)  < size(A, dim=2)) V_H = A_
  end subroutine svd_divide_and_conquer_complex_dp


  !> See [[svd_divide_and_conquer_complex_dp]]. Use this interface, if \(\mathbf{A}\) can be mutated.
  !>
  !> This routine acts on the arrays as expected by [[zgesdd]].
  subroutine zgesdd_wrapper(jobz, A, sigma, U, V_H, tol)
    !> Specifies options for computing all or part of the matrix U (see zgesdd):
    !>
    !> - `'A'`: All \( m \) left and all \( n \) right singular values are calculated.
    !>   The left singular vectors are stored column-wise in **U**.
    !>   The complex-conjugate of the right singular vectors are stored row-wise in **V_H**.
    !>
    !> - `'S'`: The first \( \min(m, n) \) clumns of **U** and rows of **V_H** are calculated
    !>
    !> - `'O'`:
    !>
    !>     - If \( m \ge n \):
    !>
    !>         - The first \( n \) columns of **A** are overwritten with  the first \( n \)
    !>           left singular vectors.
    !>
    !>         - **U** is not referenced.
    !>
    !>         - **V_H** contains the complex-conjugate of all right singular vectors, stored row-wise.
    !>
    !>     - If \( m < n \):
    !>
    !>         - The first \( m \) rows of **A** are overwritten with the complex conjugate of
    !>           the first \( m \) right singular vectors.
    !>
    !>         - **U** contains all left singular vectors, stored column-wise.
    !>
    !>         - **V_H** is not referenced.
    !>
    !> - `'N'`:  Neither left nor right singular vectors are calculated.
    character :: jobz
    !> On entry, contains the matrix \( \mathbf{A} \in \mathbb{C}^{m \times n}\).
    !>
    !> On exit:
    !>
    !> - If `jobz == 'O'`:
    !>
    !>     - If \( m \ge n \), the first \( n \) columns of **A** are overwritten with  the first \( n \)
    !>       left singular vectors.
    !>
    !>     - If \( m < n \), the first \( m \) rows of **A** are overwritten with the complex-conjugate of
    !>       the first \( m \) right singular vectors.
    !>
    !> - If `jobz \= 'O'`, the contents of A are destroyed.
    complex(dp), intent(inout), contiguous :: A(:, :)
    !> Singular values (diagonal elements of \( \Sigma \)), stored as array such that \(\sigma_i \ge \sigma_{i+1}\)
    real(dp), intent(out), contiguous :: sigma(:)
    !> Contains (part of) the left singular vectors \( \mathbf{U} \in \mathbb{C}^{m \times m} \).
    !> The vectors are stored column-wise in the array **U**.
    !>
    !> - If `jobz == 'A'`, **U** contains all left singular vectors, stored column-wise.
    !>
    !> - If `jobz == 'S'`, **U** contains the first \( \min(m, n) \) left singular vectors, stored column-wise.
    !>
    !> - If `jobz == 'O'`:
    !>
    !>     - If \( m \ge n \), **U** is not referenced.
    !>
    !>     - If \( m < n \), **U** contains all left singular vectors, stored column-wise.
    !>
    !> - If `jobz == 'N'`, **U** is not referenced.
    complex(dp), intent(out), contiguous :: U(:,:)
    !> Contains (part of) the complex-transpose of the right singular vectors \( \mathbf{V}^\dagger \in \mathbb{C}^{n \times n} \).
    !> The complex-conjugate of the vectors are stored row-wise in the array **V_H**.
    !>
    !> - If `jobz == 'A'`, **V_H** contains all left singular vectors, stored column-wise.
    !>
    !> - If `jobz == 'S'`, **V_H** contains the first \( \min(m, n) \) left singular vectors, stored column-wise.
    !>
    !> - If `jobz == 'O'`:
    !>
    !>     - If \( m \ge n \), **V_H** contains the complex-conjugate of all right singular vectors, stored row-wise.
    !>
    !>     - If \( m < n \), **V_H** is not referenced.
    !>
    !> - If `jobz == 'N'`, **V_H** is not referenced.
    complex(dp), intent(out), contiguous :: V_H(:,:)
    !> Tolerance \( \frac {\min(m, n)} {\max(m, n)} < \text{tol} \) such that 
    !> \( {\min(m, n)}  \ll {\max(m, n)} \) can be assumed
    real(dp), intent(in), optional :: tol 
    
    integer :: m, n, l, k, info, lwork, rwork_size
    real(dp) :: tol_
    integer, allocatable :: iwork(:)
    real(dp), allocatable :: rwork(:) 
    complex(dp), allocatable :: work(:)

    tol_ = default_tol_narrow_matrix
    if (present(tol)) tol_ = tol

    m = size(A, dim=1) ! rows
    n = size(A, dim=2) ! columns
    k = min(m, n)
    l = max(m, n)

    call assert_array_shape(jobz, m, n, k, size(sigma), shape(U), shape(V_H))

    allocate(work(1))
    allocate(iwork(8 * k))
    
    if (jobz == 'N') then
      rwork_size = 5*k
    else if (k / l < tol_) then
      rwork_size = 5*k**2 + 5*k 
    else
      rwork_size = max(5*k**2 + 5*k, 2*k*l + 2*k**2 + k)
    end if 
    allocate(rwork(rwork_size))

    call zgesdd(jobz, m, n, A, m, sigma, U, m, V_H, n, work, -1, rwork, iwork, info)
    call terminate_if_false(info == 0, &
                           'zgesdd failed for searching optiomal work space.')

    lwork = nint(real(work(1)))
    deallocate(work); allocate(work(lwork))

    call zgesdd(jobz, m, n, A, m, sigma, U, m, V_H, n, work, lwork, rwork, iwork, info)
    
    call terminate_if_false(info == 0, &
                           'zgesdd failed for calculating the SVD of A.')
  end subroutine zgesdd_wrapper

! Utils

  !> Setup the variable jobz for the *gesdd calls from the array sizes.
  function setup_jobz(shape_A, shape_U, shape_V_T) result(jobz)
    !> Shape of the input matrix A
    integer, intent(in) :: shape_A(2)
    !> Shape of the output matrices \(\mathbf{U}\) and  \(\mathbf{V}^\top\)
    integer, intent(in) :: shape_U(2), shape_V_T(2)

    character(len=1) jobz

    integer :: m, n, k

    m = shape_A(1)
    n = shape_A(2)
    k = min(m, n)

    if (all(shape_U == [m, m]) .and. all(shape_V_T == [n, n])) then
      jobz = 'A'   

    else if (m >= n .and. all(shape_U == shape_A) .and. all(shape_V_T == [n, n])  ) then
      jobz = 'O'

    else if (m < n .and. all(shape_U == [m, m])  .and. all(shape_V_T == shape_A) ) then
      jobz = 'O'
    
    else ! Error code
      jobz = 'X'
    end if

  end function setup_jobz

  !> Assert if the shape of the output arrays for *gesdd are valid.
  subroutine assert_array_shape(jobz, m, n, k, size_sigma, shape_U, shape_V_T)
    !> Job that shall be calculated
    character(len=1), intent(in) :: jobz
    !> Number of rows and columns of \( \mathbf{A} \) and their minimum
    integer, intent(in) :: m, n, k
    !> Size of the array containing the diagonal elements of \( \Sigma \)
    integer, intent(in) :: size_sigma
    !> Shape of \(\mathbf{U}\)
    integer, intent(in) :: shape_U(2)
    !> shape of  \(\mathbf{V}^\top\) or \(\mathbf{V}^\dagger\) respectively
    integer, intent(in) :: shape_V_T(2)

! This preprocessor usage is deliberate, to prevent if statements being evaluated in production code.
#ifdef USE_ASSERT
    call assert(jobz /= 'X', 'Sizes of arrays do not fit.')

    call assert(size_sigma == k, 'The size of sigma must be equal to min(m, n), &
                                  when A is a m by n matrix')

    if (any(jobz == ['N', 'N'])) then
      continue

    else if (any(jobz == ['A', 'A'])) then
      call assert(all(shape_U == [m, m]), &
                 'If jobz = "A", U must be a quadratic matrix with size m by m, &
                 when A is a m by n matrix')
      call assert(all(shape_V_T == [n, n]), &
                 'If jobz = "A", V_T must be a quadratic matrix with size n by n, &
                 when A is a m by n matrix')

    else if (any(jobz == ['S', 'S'])) then
      call assert(all(shape_U == [m, k]), &
                 'If jobz = "S", U must be a matrix with size m by min(m, n), &
                 when A is a m by n matrix')
      call assert(all(shape_V_T == [m, k]), &
                 'If jobz = "S", V_T must be a matrix with size min(m, n) by n, &
                 when A is a m by n matrix')

    else if (any(jobz == ['O', 'O']) .and. m >= n) then
      call assert(all(shape_U == [m, n]), &
                 'If jobz = "O" and m >= n, U must have the same shape a A.')
      call assert(all(shape_V_T == [n, n]), &
                 'If jobz = "O" and m >= n, V_T must be a quadratic matrix with size n by n, &
                 when A is a m by n matrix')

    else if (any(jobz == ['O', 'O']) .and. m < n) then
      call assert(all(shape_U == [m, m]), &
                 'If jobz = "O" and m < n, U must be a quadratic matrix with size m by m, &
                 when A is a m by n matrix')
      call assert(all(shape_V_T == [m, n]), &
                 'If jobz = "O" and m < n, V_T must have the same shape a A.')

    else
      call assert(.false., 'jobz must be one of "A", "a", "S", "s", "O", "o", "N", "n".')
    end if
#endif
  end subroutine assert_array_shape
  
end module singular_value_decomposition
