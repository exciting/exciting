!> Interfaces for singular value decomposition (svd) and matrix rank calculation, 
!> relying on svd. The interface for svd combines wrappers for the LAPACK routines
!> DGESDD, ZGESDD
module singular_value_decomposition
  use precision, only: dp
  use constants, only: zzero, zone
  use asserts, only: assert
  use modmpi, only: mpiglobal, terminate_if_false
  use math_utils, only: is_square, all_zero
  use lapack_f95_interfaces, only: dgesdd, zgesdd

  implicit none 
  
  private
  public :: svd_divide_conquer, matrix_rank

  !> Tolerance for **min(n_rows, n_cols) / max(n_rows, n_cols)** 
  real(dp), parameter :: default_tol_narrow_matrix = 1e-5_dp
  !> Tolerance for defining \(\sigma(i) = 0 \)
  real(dp), parameter :: default_tol_sigma = 1e-10_dp

  !> Calculate the singular value decomposition (SVD) of a matrix \( \mathbf{A} \):
  !> \[
  !>    \mathbf{A} = \mathbf{U} \cdot \sigma \cdot \mathbf{V}^T
  !> \]
  !> where \( \sigma \) is an m-by-n matrix which is zero except for its
  !> \( \min(m,n) \) diagonal elements, \( U \) is an \( m \)-by-\( m \) orthogonal matrix, and
  !> \( V \) is an \( n \)-by-\( n \) orthogonal matrix.  The diagonal elements of \( \sigma \)
  !> are the singular values of \( A \); they are real and non-negative, and
  !> are returned in descending order.  The first \( \min(m,n) \) columns of
  !> \( U \) and \( V \) are the left and right singular vectors of \( A \).
  !>
  !> Note that the routine returns \( \mathbf{V}^T \), not \( V \).
  interface svd_divide_conquer
    module procedure svd_divide_and_conquer_real_dp_immutable, &
                     svd_divide_and_conquer_real_dp_mutable, &
                     svd_divide_and_conquer_complex_dp_immutable, &
                     svd_divide_and_conquer_complex_dp_mutable
  end interface svd_divide_conquer


  !> Determine the rank of a matrix \( \mathbf{A} \) by calculating the SVD and 
  !> counting the non-zero singular values.
  interface matrix_rank 
    module procedure matrix_rank_real_dp, matrix_rank_complex_dp
  end interface matrix_rank

  contains 

! singular value decomposition with divide and conquer algorithm

  !> Calculate the singular value decomposition (SVD) of a matrix \( A \in \mathcal{R}^{n \times n} \).
  !> The SVD is written as
  !> \[
  !>    A = \mathbf{U} \cdot \sigma \cdot \mathbf{V}^T
  !> \]
  !> where \( \sigma \) is an m-by-n matrix which is zero except for its
  !> \( \min(m,n) \) diagonal elements, \( U \) is an \( m \)-by-\( m \) orthogonal matrix, and
  !> \( V \) is an \( n \)-by-\( n \) orthogonal matrix.  The diagonal elements of \( \sigma \)
  !> are the singular values of \( A \); they are real and non-negative, and
  !> are returned in descending order.  The first \( \min(m,n) \) columns of
  !> \( U \) and \( V \) are the left and right singular vectors of \( A \).
  !>
  !> Note that the routine returns \( \mathbf{V}^T \), not \( V \).
  !>
  !> The divide and conquer algorithm makes very mild assumptions about
  !> floating point arithmetic. It will work on machines with a guard
  !> digit in add/subtract, or on those binary machines without guard
  !> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
  !> Cray-2. It could conceivably fail on hexadecimal or decimal machines
  !> without guard digits, but we know of none.
  !> Use this interface, if \(\mathbf{A}\) should be unchanged.
  !> (This documation is from netlib.org)
  subroutine svd_divide_and_conquer_real_dp_immutable(A, sigma, U, V_T)
    !> Input matrix \( A \)
    real(dp), intent(in), contiguous :: A(:, :)
    !> Singular values \( \sigma \)
    real(dp), intent(out), contiguous :: sigma(:)
    !> Left singular vectors \( \mathbf{U} \)
    !> <li> If not given, only \( \sigma \) is calculated.
    !> <li> If \( \mathbf{U} \in \mathcal{C}^{m \times m} \), 
    !> all left singular vectors are calculated
    !> <li> If \( \mathbf{U} \in \mathcal{C}^{m \times n} \) or 
    !> \( \in \mathcal{C}^{n \times m} \) and \(m >= n \), 
    !> the first 'n' columns of \( \mathbf{U} \) are calculated.
    real(dp), intent(out), contiguous, optional :: U(:,:)
    !> Adjungate of right singular vectors \( \mathbf{V}^T \)
    !> Complex-transpose of the right singular vectors \( \mathbf{V}^T \)
    !> <li> If not given, only \( \sigma \) is calculated.
    !> <li> If \( \mathbf{V}^T \in \mathcal{C}^{n \times n} \), 
    !> all right singular vectors are calculated
    !> <li> If \( \mathbf{V}^T \in \mathcal{C}^{m \times n} \) and \(m >= n \), 
    !> only the first 'n' rows of \( \mathbf{V}^T \) are calculated.
    real(dp), intent(out), contiguous, optional :: V_T(:,:)
        
    character(1) :: jobz
    real(dp) :: U_(1, 1), V_T_(1, 1)
    real(dp), allocatable :: A_(:, :)
      
    ! dgesdd destroys the input
    A_ = A

    if(present(U) .and. present(V_T)) then
      jobz = setup_jobz(shape(A), shape(U), shape(V_T))
      call svd_divide_and_conquer_real_dp_mutable(jobz, A_, sigma, U, V_T)
    
    else if(.not. (present(U) .or. present(V_T))) then
      jobz = 'N'
      call svd_divide_and_conquer_real_dp_mutable(jobz, A_, sigma, U_, V_T_)
    
    else 
      call terminate_if_false(.false., 'Either both, U and V_T, must be present or none of both.')
    end if

    if (jobz == 'O' .and. size(A, dim=1) >= size(A, dim=2)) U = A_ 
    if (jobz == 'O' .and. size(A, dim=1)  < size(A, dim=2)) V_T = A_ 
  end subroutine svd_divide_and_conquer_real_dp_immutable


  !> See [[svd_divide_and_conquer_real_dp_immutable]]. Use this interface if \(\mathbf{A}\) can be mutated.
  subroutine svd_divide_and_conquer_real_dp_mutable(jobz, A, sigma, U, V_T)
    !> Specifies options for computing all or part of the matrix U:
    !> <li>      = 'A':  all m columns of \(\mathbf{U}\) and all n rows of  \(\mathbf{V}^T\)  are
    !>               returned in the arrays \(\mathbf{U}\) and  \(\mathbf{V}^T\) ;
    !> <li>      = 'S':  the first \(\min(m, n)\) columns of \(\mathbf{U}\) and the first
    !>               \(\min(m, n)\) rows of  \(\mathbf{V}^T\)  are returned in the arrays U
    !>               and  \(\mathbf{V}^T\) ;
    !> <li>      = 'O':  If m >= n, the first n columns of \(\mathbf{U}\) are overwritten
    !>               in the array A and all rows of  \(\mathbf{V}^T\)  are returned in
    !>               the array  \(\mathbf{V}^T\) ;
    !>               otherwise, all columns of \(\mathbf{U}\) are returned in the
    !>               array \(\mathbf{U}\) and the first m rows of  \(\mathbf{V}^T\)  are overwritten
    !>               in the array A;
    !> <li>      = 'N':  no columns of \(\mathbf{U}\) or rows of  \(\mathbf{V}^T\)  are computed.
    character :: jobz
    !> On entry, input matrix \( A \in \mathcal{R}^{m \times n} \).
    !>      On exit,
    !> <li>     if jobz = 'O',  A is overwritten with the first 'n' columns
    !>                      of \(\mathbf{U}\) (the left singular vectors, stored
    !>                      columnwise) if \(m >= n\);
    !>                      A is overwritten with the first 'm' rows
    !>                      of V_H (the right singular vectors, stored
    !>                      rowwise) otherwise.
    !> <li>     if jobz \= 'O', the contents of A are destroyed.
    real(dp), intent(inout), contiguous :: A(:, :)
    !> Singular values \( \sigma \), sorted such that \(\sigma_i = \sigma_{i+1}\)
    real(dp), intent(out), contiguous :: sigma(:)
    !> \(\mathbf{U}\) is a real \( {m \times n_U} \) array. 
    !>       \(n_U = m\) if jobz = 'A' or jobz = 'O' and  \(m >= n\);
    !>       \(n_U = \min(m, n) if jobz = 'S'.
    !> <li>      If jobz = 'A' or jobz = 'O' and  \(m >= n\), \(\mathbf{U}\) contains the unitary matrix 
    !>                      \( \mathbf{U}  \in \mathcal{R}^{m \times n} \);
    !> <li>      if jobz = 'S', \(\mathbf{U}\) contains the first \(\min(m, n)\) columns of \( \mathbf{U} \)
    !>                      (the left singular vectors, stored columnwise);
    !> <li>      if jobz = 'O' and m >= n, or jobz = 'N', \(\mathbf{U}\) is not referenced.
    real(dp), intent(out), contiguous :: U(:,:)
    !>  \(\mathbf{V}^T\)  is a real \( {n \times n_V} \) array.
    !> <li>      If jobz = 'A' or jobz = 'O' and \(m >= n\),  \(\mathbf{V}^T\)  contains the the unitary matrix
    !>                      \(\mathbf{V}^T \in \mathcal{R}^{n \times n}\);
    !> <li>      if jobz = 'S',  \(\mathbf{V}^T\)  contains the first \(\min(m, n)\) rows of
    !>                      \(\mathbf{V}^T\) (the right singular vectors, stored rowwise);
    !> <li>      if jobz = 'O' and m < n, or jobz = 'N',  \(\mathbf{V}^T\)  is not referenced.
    real(dp), intent(out), contiguous :: V_T(:,:)
    
    integer :: n_rows, n_cols, k, info, lwork
    integer, allocatable :: iwork(:)
    real(dp), allocatable :: work(:)


    n_rows = size(A, dim=1)
    n_cols = size(A, dim=2)
    k = min(n_rows, n_cols)

    call assert_array_shape(jobz, n_rows, n_cols, k, size(sigma), shape(U), shape(V_T))
    
    lwork = -1
    allocate(work(1))
    allocate(iwork(8 * k))
      
    call dgesdd(jobz, n_rows, n_cols, A, n_rows, sigma, U, n_rows, V_T, n_cols, work, lwork, iwork, info)
    call terminate_if_false(info == 0, &
                           'dgesdd failed for searching optiomal work space.')

    lwork = idnint(work(1))
    deallocate(work); allocate(work(lwork))
    
    call dgesdd(jobz, n_rows, n_cols, A, n_rows, sigma, U, n_rows, V_T, n_cols, work, lwork, iwork, info)

    call terminate_if_false(info == 0, &
                           'dgesdd failed for calculating the SVD of A.')
  end subroutine svd_divide_and_conquer_real_dp_mutable


  !> Calculate the singular value decomposition (SVD) of a matrix 
  !> \( \mathbf{A} \in \mathcal{C}^{n \times n} \).
  !> The SVD is written as
  !> \[
  !>    \mathbf{A} = \mathbf{U} \cdot \sigma \cdot V^\dagger
  !> \]
  !> where \( \sigma \) is an m-by-n matrix which is zero except for its
  !> \( \min(m,n) \) diagonal elements, \( \mathbf{U} \) is an \( m \)-by-\( m \) orthogonal matrix, and
  !> \( V \) is an \( n \)-by-\( n \) orthogonal matrix.  The diagonal elements of \( \sigma \)
  !> are the singular values of \( A \); they are real and non-negative, and
  !> are returned in descending order.  The first \( \min(m,n) \) columns of
  !> \( \mathbf{U} \) and \( \mathbf{V} \) are the left and right singular vectors of \( A \).
  !>
  !> Note that the routine returns \( \mathbf{V}^\dagger \), not \( V \).
  !>
  !> The divide and conquer algorithm makes very mild assumptions about
  !> floating point arithmetic. It will work on machines with a guard
  !> digit in add/subtract, or on those binary machines without guard
  !> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
  !> Cray-2. It could conceivably fail on hexadecimal or decimal machines
  !> without guard digits, but we know of none.
  !> Use this interface, if \(\mathbf{A}\) should be unchanged.
  !> (This documation is from netlib.org)
  subroutine svd_divide_and_conquer_complex_dp_immutable(A, sigma, U, V_C)
    !> Input matrix \( \mathbf{A} \in \mathcal{C}^{m \times n} \)
    complex(dp), intent(in), contiguous :: A(:,:)
    !> Singular values \( \sigma \)
    real(dp), intent(out), contiguous :: sigma(:)
    !> Left singular vectors \( \mathbf{U} \)
    !> <li> If not given, only \( \sigma \) is calculated.
    !> <li> If \( \mathbf{U} \in \mathcal{C}^{m \times m} \), 
    !> all left singular vectors are calculated
    !> <li> If \( \mathbf{U} \in \mathcal{C}^{m \times n} \) or 
    !> \( \in \mathcal{C}^{n \times m} \) and \(m >= n \), 
    !> the first 'n' columns of \( \mathbf{U} \) are calculated.
    complex(dp), intent(out), contiguous, optional :: U(:,:)
    !> Adjungate of right singular vectors \( \mathbf{V}^\dagger \)
    !> Complex-transpose of the right singular vectors \( \mathbf{V}^\dagger \)
    !> <li> If not given, only \( \sigma \) is calculated.
    !> <li> If \( \mathbf{V}^\dagger \in \mathcal{C}^{n \times n} \), 
    !> all right singular vectors are calculated
    !> <li> If \( \mathbf{V}^\dagger \in \mathcal{C}^{m \times n} \) and \(m >= n \), 
    !> only the first 'n' rows of \( \mathbf{V}^\dagger \) are calculated.
    complex(dp), intent(out), contiguous, optional :: V_C(:,:)
        
    character(1) :: jobz
    complex(dp) :: U_(1, 1), V_C_(1, 1)
    complex(dp), allocatable :: A_(:, :)
      
    ! zgesdd destroys the input
    A_ = A

    if(present(U) .and. present(V_C)) then
      jobz = setup_jobz(shape(A), shape(U), shape(V_C))
      call svd_divide_and_conquer_complex_dp_mutable(jobz, A_, sigma, U, V_C)
    else
      jobz = 'N'
      call svd_divide_and_conquer_complex_dp_mutable(jobz, A_, sigma, U_, V_C_)
    end if

    if (jobz == 'O' .and. size(A, dim=1) >= size(A, dim=2)) U = A_ 
    if (jobz == 'O' .and. size(A, dim=1)  < size(A, dim=2)) V_C = A_ 
  end subroutine svd_divide_and_conquer_complex_dp_immutable


  !> See [[svd_divide_and_conquer_complex_dp_immutable]]. Use this interface, if \(\mathbf{A}\) can be mutated.
  subroutine svd_divide_and_conquer_complex_dp_mutable(jobz, A, sigma, U, V_C, tol)
    !> Specifies options for computing all or part of the matrix U:
    !> <li>      = 'A':  all m columns of \(\mathbf{U}\) and all n rows of  \(\mathbf{V}^\dagger\)  are
    !>               returned in the arrays \(\mathbf{U}\) and  \(\mathbf{V}^\dagger\) ;
    !> <li>      = 'S':  the first \(\min(m, n)\) columns of \(\mathbf{U}\) and the first
    !>               \(\min(m, n)\) rows of  \(\mathbf{V}^\dagger\)  are returned in the arrays U
    !>               and  \(\mathbf{V}^\dagger\) ;
    !> <li>      = 'O':  If m >= n, the first n columns of \(\mathbf{U}\) are overwritten
    !>               in the array A and all rows of  \(\mathbf{V}^\dagger\)  are returned in
    !>               the array  \(\mathbf{V}^\dagger\) ;
    !>               otherwise, all columns of \(\mathbf{U}\) are returned in the
    !>               array \(\mathbf{U}\) and the first m rows of  \(\mathbf{V}^\dagger\)  are overwritten
    !>               in the array A;
    !> <li>      = 'N':  no columns of \(\mathbf{U}\) or rows of  \(\mathbf{V}^\dagger\)  are computed.
    character :: jobz
    !> On entry, input matrix \( A \in \mathcal{R}^{m \times n}\).
    !>      On exit,
    !> <li>     if jobz = 'O',  A is overwritten with the first 'n' columns
    !>                      of \(\mathbf{U}\) (the left singular vectors, stored
    !>                      columnwise) if \(m >= n\);
    !>                      A is overwritten with the first 'm' rows
    !>                      of V_H (the right singular vectors, stored
    !>                      rowwise) otherwise.
    !> <li>     if jobz \= 'O', the contents of A are destroyed.
    complex(dp), intent(inout) :: A(:, :)
    !> Singular values \( \sigma \), sorted such that \(\sigma_i = \sigma_{i+1}\)
    real(dp), intent(out) :: sigma(:)
    !> \(\mathbf{U}\) is a complex \( {m \times n_U} \) array. 
    !>       \(n_U = m\) if jobz = 'A' or jobz = 'O' and  \(m >= n\);
    !>       \(n_U = \min(m, n) if jobz = 'S'.
    !> <li>      If jobz = 'A' or jobz = 'O' and  \(m >= n\), \(\mathbf{U}\) contains the unitary matrix 
    !>                      \( \mathbf{U}  \in \mathcal{R}^{m \times n} \);
    !> <li>      if jobz = 'S', \(\mathbf{U}\) contains the first \(\min(m, n)\) columns of \( \mathbf{U} \)
    !>                      (the left singular vectors, stored columnwise);
    !> <li>      if jobz = 'O' and m >= n, or jobz = 'N', \(\mathbf{U}\) is not referenced.
    complex(dp), intent(out) :: U(:,:)
    !>  \(\mathbf{V}^\dagger\)  is a complex \( {n \times n_V} \) array.
    !> <li>      If jobz = 'A' or jobz = 'O' and \(m >= n\),  \(\mathbf{V}^\dagger\)  contains the the unitary matrix
    !>                      \(\mathbf{V}^T \in \mathcal{R}^{n \times n}\);
    !> <li>      if jobz = 'S',  \(\mathbf{V}^\dagger\)  contains the first \(\min(m, n)\) rows of
    !>                      \(\mathbf{V}^T\) (the right singular vectors, stored rowwise);
    !> <li>      if jobz = 'O' and m < n, or jobz = 'N',  \(\mathbf{V}^\dagger\)  is not referenced.
    complex(dp), intent(out) :: V_C(:,:)
    !> Tolerance \( \frac {\min(m, n)} {\max(m, n)} < \text{tol} \) such that 
    !> \( {\min(m, n)}  \ll {\max(m, n)} \) can be assumed
    real(dp), intent(in), optional :: tol 
    
    integer :: n_rows, n_cols, l, k, info, lwork, rwork_size
    real(dp) :: tol_
    integer, allocatable :: iwork(:)
    real(dp), allocatable :: rwork(:) 
    complex(dp), allocatable :: work(:)

    tol_ = default_tol_narrow_matrix
    if (present(tol)) tol_ = tol

    n_rows = size(A, dim=1) ! rows
    n_cols = size(A, dim=2) ! columns
    k = min(n_rows, n_cols)
    l = max(n_rows, n_cols)

    call assert_array_shape(jobz, n_rows, n_cols, k, size(sigma), shape(U), shape(V_C))

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

    call zgesdd(jobz, n_rows, n_cols, A, n_rows, sigma, U, n_rows, V_C, n_cols, work, -1, rwork, iwork, info)
    call terminate_if_false(info == 0, &
                           'zgesdd failed for searching optiomal work space.')

    lwork = nint(real(work(1)))
    deallocate(work); allocate(work(lwork))

    call zgesdd(jobz, n_rows, n_cols, A, n_rows, sigma, U, n_rows, V_C, n_cols, work, lwork, rwork, iwork, info)
    
    call terminate_if_false(info == 0, &
                           'zgesdd failed for calculating the SVD of A.')
  end subroutine svd_divide_and_conquer_complex_dp_mutable


! matrix_rank

  !> Calculate the rank of a matrix \( A \) by determining
  !> the singular value decomposition. The rank is 
  !> given by the number of non-zero singular values
  function matrix_rank_real_dp(A, tol) result(rank)
    !> Input matrix \( A \)
    real(dp), intent(in), contiguous :: A(:, :)
    !> Tolerance for defining \(\sigma(i) = 0\)
    real(dp), intent(in), optional :: tol 
  
    integer :: rank

    real(dp) :: tol_
    real(dp), allocatable :: sigma(:)
    integer :: i

    tol_ = default_tol_sigma
    if(present(tol)) tol_ = tol

    allocate(sigma(min(size(A, dim=1), size(A, dim=2))))
    call svd_divide_conquer(A, sigma)

    rank = 0
    do i=1, size(sigma)
      if (.not. all_zero(sigma(i))) rank = rank + 1
    end do
  end function matrix_rank_real_dp


  !> Calculate the rank of a matrix \( A \) by determining
  !> the singular value decomposition. The rank is 
  !> given by the number of non-zero singular values.
  function matrix_rank_complex_dp(A, tol) result(rank)
    !> Input matrix \( A \)
    complex(dp), intent(in), contiguous :: A(:, :)
    !> Tolerance for defining \(\sigma(i) = 0\)
    real(dp), intent(in), optional :: tol 
  
    integer :: rank

    real(dp) :: tol_
    real(dp), allocatable :: sigma(:)
    integer :: i

    tol_ = default_tol_sigma
    if(present(tol)) tol_ = tol

    allocate(sigma(min(size(A, dim=1), size(A, dim=2))))
    call svd_divide_conquer(A, sigma)

    rank = 0
    do i=1, size(sigma)
      if (.not. all_zero(sigma(i))) rank = rank + 1
    end do
  end function matrix_rank_complex_dp


! Utils

  !> Setup the variable jobz for the *gesdd calls from the array sizes.
  function setup_jobz(shape_A, shape_U, shape_V_T) result(jobz)
    !> Shape of the input matrix A
    integer, intent(in) :: shape_A(2)
    !> Shape of the output matrices \(\mathbf{U}\) and  \(\mathbf{V}^T\)  
    integer, intent(in) :: shape_U(2), shape_V_T(2)

    character(len=1) jobz

    integer :: n_rows, n_cols, k 

    n_rows = shape_A(1)
    n_cols = shape_A(2)
    k = min(n_rows, n_cols)

    if (all(shape_U == [n_rows, n_rows]) .and. all(shape_V_T == [n_cols, n_cols])) then 
      jobz = 'A'   

    else if (n_rows >= n_cols .and. all(shape_U == shape_A) .and. all(shape_V_T == [n_cols, n_cols])  ) then
      jobz = 'O'

    else if (n_rows < n_cols .and. all(shape_U == [n_rows, n_rows])  .and. all(shape_V_T == shape_A) ) then
      jobz = 'O'
    
    else ! Error code
      jobz = 'X'
    end if

  end function setup_jobz


  !> Assert if the shape of the output arrays for *gesdd are valid.
  subroutine assert_array_shape(jobz, n_rows, n_cols, k, size_sigma, shape_U, shape_V_T)
    !> Job that shall be calculated
    character(len=1), intent(in) :: jobz
    !> Number of rows and columns of A and the minimum
    integer, intent(in) :: n_rows, n_cols, k 
    !> Size of sigma, shape of \(\mathbf{U}\) and shape of  \(\mathbf{V}^T\) 
    integer, intent(in) :: size_sigma, shape_U(2), shape_V_T(2)

#ifdef USE_ASSERT
    call assert(jobz /= 'X', 'Sizes of arrays do not fit.')

    call assert(size_sigma == k, 'The size of sigma must be equal to min(n_rows, n_cols), &
                                  when A is a n_rows by n_cols matrix')

    if (any(jobz == ['n', 'N'])) then
      continue

    else if (any(jobz == ['a', 'A'])) then
      call assert(all(shape_U == [n_rows, n_rows]), &
                 'If jobz = "A", U must be a quadratic matrix with size n_rows by n_rows, &
                 when A is a n_rows by n_cols matrix')
      call assert(all(shape_V_T == [n_cols, n_cols]), &
                 'If jobz = "A", V_T must be a quadratic matrix with size n_cols by n_cols, &
                 when A is a n_rows by n_cols matrix')

    else if (any(jobz == ['s', 'S'])) then 
      call assert(all(shape_U == [n_rows, k]), &
                 'If jobz = "S", U must be a matrix with size n_rows by min(n_rows, n_cols), &
                 when A is a n_rows by n_cols matrix')
      call assert(all(shape_V_T == [n_rows, k]), &
                 'If jobz = "S", V_T must be a matrix with size min(n_rows, n_cols) by n_cols, &
                 when A is a n_rows by n_cols matrix')

    else if (any(jobz == ['o', 'O']) .and. n_rows >= n_cols) then
      call assert(all(shape_U == [n_rows, n_cols]), &
                 'If jobz = "O" and n_rows >= n_cols, U must have the same shape a A.')
      call assert(all(shape_V_T == [n_cols, n_cols]), &
                 'If jobz = "O" and n_rows >= n_cols, V_T must be a quadratic matrix with size n_cols by n_cols, &
                 when A is a n_rows by n_cols matrix')

    else if (any(jobz == ['o', 'O']) .and. n_rows < n_cols) then
      call assert(all(shape_U == [n_rows, n_rows]), &
                 'If jobz = "O" and n_rows < n_cols, U must be a quadratic matrix with size n_rows by n_rows, &
                 when A is a n_rows by n_cols matrix')
      call assert(all(shape_V_T == [n_rows, n_cols]), &
                 'If jobz = "O" and n_rows < n_cols, V_T must have the same shape a A.')    

    else
      call assert(.false., 'jobz must be one of "A", "a", "S", "s", "O", "o", "N", "n".')
    end if
#endif

  end subroutine assert_array_shape
  
end module singular_value_decomposition