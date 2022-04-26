!> Interfaces for LU factorization and routines that use LU factorizations.
!> The interfaces combine wrappers for the LAPACK routines 
!> **[[dgetc2]]**, **[[dgetrf2]]**, **[[zgetc2]]**, **[[zgetrf2]]**, **[[dgetri]]**, **[[zgetri]]**
module lu_factorization
  use asserts, only: assert
  use modmpi, only: terminate_if_false
  use precision, only: dp
  use constants, only: zzero, zone
  use math_utils, only: is_square
  use grid_utils, only: mesh_1d
  use lapack_f95_interfaces, only: dgetc2, zgetc2, dgetrf2, zgetrf2, dgetri, zgetri
  use matrix_rank, only: matrix_rank_SVD

  private
  public :: lu_full_pivot, lu_row_pivot, xgetc2, xgetrf2, xgetri, pivot_to_permutation


  !> Default of `caclulate_permutation`, which specifies if the permutation of the rows and columns is
  !> calculated from the pivot vectors as returned by the LAPACK routines. If false, the pivot vector is returned.
  !> As this is an additional calculation on top of the LAPACK routine, the default is `.false.`.
  logical, parameter :: default_calculate_permutation = .false.


  !> Calculate the LU factorization with full pivoting for a quadratic matrix \( \mathbf{A} \)
  !> such that
  !> \[
  !>   \mathbf{A} = \mathbf{P} \cdot \mathbf{L} \cdot \mathbf{U},
  !> \]
  !> where \( \mathbf{P} \) is the row permutation matrix, \( \mathbf{L} \) is a lower triangular 
  !> matrix with unit diagonal elements and \( \mathbf{U} \) is upper triangular.
  !> 
  !> The interface does not change \( \mathbf{A} \) and writes the results in the arrays \( \mathbf{L} \) and 
  !> \( \mathbf{U} \) and optionally allows for calculating the permutation vectors from the pivot vectors.
  !> 
  !> For more detailed information see the documentation from **[[lu_full_pivot_real_dp]]** and 
  !> **[[lu_full_pivot_complex_dp]]**.
  interface lu_full_pivot
    module procedure lu_full_pivot_real_dp, &
                     lu_full_pivot_complex_dp
  end interface lu_full_pivot


  !> See **[[lu_full_pivot]]**. 
  !>
  !> \( \mathbf{A} \) will be overwritten with this interface and only the pivot vectors will be returned.
  !>
  !> For more detailed information see the documentation from **[[dgetc2_wrapper]]** and 
  !> **[[zgetc2_wrapper]]**.
  interface xgetc2
    module procedure dgetc2_wrapper, &
                     zgetc2_wrapper
  end interface xgetc2


  !> Calculate the LU factorization with row pivoting for a quadratic matrix \( \mathbf{A} \)
  !> such that
  !> \[
  !>   \mathbf{A} = \mathbf{P} \cdot \mathbf{L} \cdot \mathbf{U} \cdot \mathbf{Q},
  !> \]
  !> where \( \mathbf{P} \) and \( \mathbf{Q} \) are permutation matrices, \( \mathbf{L} \) is a lower triangular 
  !> matrix with unit diagonal elements and \( \mathbf{U} \) is upper triangular.
  !> 
  !> The interface does not change \( \mathbf{A} \) and writes the results in the arrays \( \mathbf{L} \) and 
  !> \( \mathbf{U} \) and optionally allows for calculating the permutation vector from the pivot vector.
  !> 
  !> For more detailed information see the documentation from **[[lu_row_pivot_real_dp]]** and 
  !> **[[lu_row_pivot_complex_dp]]**.
  interface lu_row_pivot 
    module procedure lu_row_pivot_real_dp, &
                     lu_row_pivot_complex_dp
  end interface lu_row_pivot

  !> See **[[lu_row_pivot]]**. 
  !>
  !> \( \mathbf{A} \) will be overwritten with this interface and only the pivot vector but not the permutation vector
  !> can be returned.
  !>
  !> For more detailed information see the documentation from **[[dgetrf2_wrapper]]** and
  !> **[[zgetrf2_wrapper]]**.
  interface xgetrf2
    module procedure dgetrf2_wrapper, &
                     zgetrf2_wrapper
  end interface xgetrf2

  !> For a LU factorization, calculated by **[[xgetrf2]]**,
  !> \[
  !>    \mathbf{A} = \mathbf{L} \cdot \mathbf{U},
  !> \]
  !> invert \( \mathbf{U} \) and solve the system of linear equations
  !> \[
  !>    \mathbf{A}^{-1} \cdot \mathbf{L} = \mathbf{U}^{-1}.
  !> \]
  !> This routine does not assert if the matrix \( \mathbf{A} \) is invertible.
  !>
  !> This is a wrapper for the LAPACK routine dgetri.
  interface xgetri
    module procedure dgetri_wrapper, &
                     zgetri_wrapper
  end interface xgetri

  !> Copy from a matrix \( \mathbf{A} \) the lower triangular part \( \mathbf{L} \) and the upper triangular part
  !> \( \mathbf{U} \) into two arrays \( \mathbf{L} \) \( \mathbf{U} \), where the diagonal of \( \mathbf{L} \) is
  !> known to be one.
  interface extract_LU
    module procedure :: extract_LU_real_dp, extract_LU_complex_dp
  end interface extract_LU

  contains

! LU factorization with full pivoting
  
  !> Calculate LU factorization with full pivoting for a real matrix \( \mathbf{A} \in \mathbb{R}^{N \times N} \)
  !> such that
  !> \[
  !>   \mathbf{A} = \mathbf{P} \cdot \mathbf{L} \cdot \mathbf{U} \cdot \mathbf{Q},
  !> \]
  !> where \( \mathbf{P} \) and \( \mathbf{Q} \) are permutation matrices, \( \mathbf{L} \) is a lower triangular 
  !> matrix with unit diagonal elements and \( \mathbf{U} \) is upper triangular.
  !> Use this interface, if  \( \mathbf{A} \) should be unchanged.
  !>
  !> This is a wrapper for the LAPACK routine dgetc2.
  subroutine lu_full_pivot_real_dp(A_input, L, U, P, Q, calculate_permutation)
    !> Input matrix \( \mathbf{A} \)
    real(dp), intent(in), contiguous :: A_input(:, :)
    !> Lower triangular matrix \( \mathbf{L} \)
    real(dp), intent(out), contiguous :: L(:, :)
    !> Upper triangular matrix \( \mathbf{U} \) 
    real(dp), intent(out), contiguous :: U(:, :)
    !> If **calculate_permutation** is `.true.`, the permutation vector of the rows is returned,
    !> else the row pivot vector as obtained by **[[dgetc2_wrapper]]** is returned.
    integer, intent(out), contiguous :: P(:)
    !> If **calculate_permutation** is `.true.`, the permutation vector of the columns is returned,
    !> else the column pivot vector as obtained by **[[dgetc2_wrapper]]** is returned.
    integer, intent(out), contiguous :: Q(:)
    !> Return either permutation vectors or pivot vectors.
    !> If `.true.`, the permutation vectors are returned,
    !> else the pivot vectors are returned.
    !> Default is `.false.`.
    logical, intent(in), optional :: calculate_permutation 
    
    !> Pivot index maps for rows (`row_pivot_map`) and columns (`column_pivot_map`)
    integer, allocatable :: row_pivot_map(:), column_pivot_map(:)
    real(dp), allocatable :: A(:, :)
    logical :: calculate_permutation_
        
    calculate_permutation_ = default_calculate_permutation
    if (present(calculate_permutation)) calculate_permutation_ = calculate_permutation

    A = A_input 
    call dgetc2_wrapper(A, P, Q)
    
    call extract_LU_real_dp(A, L, U)

    if (calculate_permutation_) then
      row_pivot_map = P 
      column_pivot_map = Q 
      call pivot_to_permutation(row_pivot_map, P)
      call pivot_to_permutation(column_pivot_map, Q)
    end if 
  end subroutine lu_full_pivot_real_dp


  !> See **[[lu_full_pivot_real_dp]]**. Use this interface if \( \mathbf{A} \) can be overwritten.
  !>
  !> This is a wrapper for the LAPACK routine dgetc2.
  subroutine dgetc2_wrapper(A, row_pivot_map, column_pivot_map)
    !> On entry, the matrix \( \mathbf{A} \), on exit 
    !> the factors \( \mathbf{L} \) and \( \mathbf{U} \), where the 
    !> unit diagonals of \( \mathbf{L} \) are not stored.
    real(dp), intent(inout), contiguous :: A(:, :)
    !> Pivot indices; for \( 1 \le \) `i` \(  \le N \), row  `A(i, :)`  has
    !> been interchanged with row `A(row_pivot_map(i), :)`.
    integer, intent(out), contiguous :: row_pivot_map(:)
    !> Pivot indices; for \( 1 \le \) `j` \(  \le N \), column  `A(:, j)`  has
    !> been interchanged with column `A(:, row_pivot_map(j))`.
    integer, intent(out), contiguous :: column_pivot_map(:)

    integer :: info

    call assert(is_square(A), 'A is not a square matrix.')
    call assert(matrix_rank_SVD(A) == size(A, dim=1), 'A is not a full rank matrix.')
    call assert(size(row_pivot_map) == size(A, dim=1), &
            'The size of row_pivot_map is not equal to the number of rows of A.')
    call assert(size(column_pivot_map) == size(A, dim=2), &
            'The size of column_pivot_map is not equal to the number of columns of A.')

    call dgetc2(size(A, dim=1), A, size(A, dim=1), row_pivot_map, column_pivot_map, info)

    call terminate_if_false(info >= 0, 'dgetc2 failed.')
  end subroutine dgetc2_wrapper


  !> Calculate LU factorization with full pivoting for a complex matrix \( \mathbf{A} \in \mathbb{C}^{N \times N} \)
  !> such that
  !> \[
  !>   \mathbf{A} = \mathbf{P} \cdot \mathbf{L} \cdot \mathbf{U} \cdot \mathbf{Q},
  !> \]
  !> where \( \mathbf{P} \) and \( \mathbf{Q} \) are permutation matrices, \( \mathbf{L} \) is a lower triangular 
  !> matrix with unit diagonal elements and \( \mathbf{U} \) is upper triangular.
  !> Use this interface, if  \( \mathbf{A} \) should be unchanged. 
  !>
  !> This is a wrapper for the LAPACK routine zgetc2.
  subroutine lu_full_pivot_complex_dp(A_input, L, U, P, Q, calculate_permutation)
    !> Input matrix \( \mathbf{A} \)
    complex(dp), intent(in), contiguous :: A_input(:, :)
    !> Lower triangular matrix \( \mathbf{L} \)
    complex(dp), intent(out), contiguous :: L(:, :)
    !> Upper triangular matrix \( \mathbf{U} \) 
    complex(dp), intent(out), contiguous :: U(:, :)
    !> If **calculate_permutation** is `.true.`, the permutation vector of the rows is returned,
    !> else the row pivot vector as obtained by **[[zgetc2_wrapper]]** is returned.
    integer, intent(out), contiguous :: P(:)
    !> If **calculate_permutation** is `.true.`, the permutation vector of the columns is returned
    !> else the column pivot vector as obtained by **[[zgetc2_wrapper]]** is returned.
    integer, intent(out), contiguous :: Q(:)
    !> Return either permutation vectors or pivot vectors.
    !> If `.true.`, the permutation vectors are returned,
    !> else the pivot vectors are returned.
    !> Default is `.false.`.
    logical, intent(in), optional :: calculate_permutation 
    
    !> Pivot index maps for rows (`row_pivot_map`) and columns (`column_pivot_map`)
    integer, allocatable :: row_pivot_map(:), column_pivot_map(:)
    complex(dp), allocatable :: A(:,:)
    logical :: calculate_permutation_
        
    calculate_permutation_ = default_calculate_permutation
    if (present(calculate_permutation)) calculate_permutation_ = calculate_permutation

    A = A_input 
    call zgetc2_wrapper(A, P, Q)
    
    call extract_LU_complex_dp(A, L, U)

    if (calculate_permutation_) then
      row_pivot_map = P 
      column_pivot_map = Q 
      call pivot_to_permutation(row_pivot_map, P)
      call pivot_to_permutation(column_pivot_map, Q)
    end if 
  end subroutine lu_full_pivot_complex_dp

   !> See **[[lu_full_pivot_complex_dp]]**. Use this interface if \( \mathbf{A} \) can be overwritten.
  !>
  !> This is a wrapper for the LAPACK routine dgetc2.
  subroutine zgetc2_wrapper(A, row_pivot_map, column_pivot_map)
    !> On entry, the matrix \( \mathbf{A} \), on exit 
    !> the factors \( \mathbf{L} \) and \( \mathbf{U} \), where the 
    !> unit diagonals of \( \mathbf{L} \) are not stored.
    complex(dp), intent(inout), contiguous :: A(:, :)
    !> Pivot indices; for \( 1 \le \) `i` \(  \le N \), row  `A(i, :)`  has
    !> been interchanged with row `A(row_pivot_map(i), :)`.
    integer, intent(out), contiguous :: row_pivot_map(:)
    !> Pivot indices; for \( 1 \le \) `j` \(  \le N \), column  `A(:, j)`  has
    !> been interchanged with column `A(:, row_pivot_map(j))`.
    integer, intent(out), contiguous :: column_pivot_map(:)

    integer :: info

    call assert(is_square(A), 'A is not a square matrix.')
    call assert(matrix_rank_SVD(A) == size(A, dim=1), 'A is not a full rank matrix.')
    call assert(size(row_pivot_map) == size(A, dim=1), &
            'The size of row_pivot_map is not equal to the number of rows of A.')
    call assert(size(column_pivot_map) == size(A, dim=2), &
            'The size of column_pivot_map is not equal to the number of columns of A.')
    
    call zgetc2(size(A, dim=1), A, size(A, dim=1), row_pivot_map, column_pivot_map, info)

    call terminate_if_false(info >= 0, 'zgetc2 failed.')
  end subroutine zgetc2_wrapper

! LU factorization with row pivoting

  !> Calculate the LU factorization with partial row pivoting for a real matrix
  !> \( \mathbf{A} \in \mathbb{R}^{M \times N} \) such that
  !> \[
  !>   \mathbf{A} = \mathbf{P} \cdot \mathbf{L} \cdot \mathbf{U},
  !> \]
  !> where \( \mathbf{P} \) is the permutation matrix, \( \mathbf{L} \) is a lower triangular 
  !> matrix with unit diagonal elements and \( \mathbf{U} \) is upper triangular.
  !> Use this interface, if  \( \mathbf{A} \) should be unchanged. 
  !>
  !> This is a wrapper for the LAPACK routine dgetrf2.
  subroutine lu_row_pivot_real_dp(A_input, L, U, P, calculate_permutation)
    !> Input matrix \( \mathbf{A} \in \mathbb{R}^{m \times n} \)
    real(dp), intent(in), contiguous :: A_input(:, :)
    !> Lower triangular matrix  \( \mathbf{L} \) (or lower trapezoidal if \(m > n\)).
    real(dp), intent(out), contiguous :: L(:, :)
    !> Upper triangular matrix \( \mathbf{U} \) (or upper trapezoidal if \(m < n\)).
    real(dp), intent(out), contiguous :: U(:, :)
    !> If **calculate_permutation** is `.true.`, the permutation vector of the rows is returned,
    !> else the row pivot vector as obtained by **[[dgetrf2]]** is returned.
    integer, intent(out), contiguous :: P(:)
    !> Return either permutation vectors or pivot vectors.
    !> If `.true.`, the permutation vectors are returned,
    !> else the pivot vectors are returned. If the pivot vector is returned, consider only the
    !> first min'(m, n)' elements.
    !> Default is `.false.`.
    logical, intent(in), optional :: calculate_permutation 

    integer :: k
    integer, allocatable :: row_pivot_map(:)
    real(dp), allocatable :: A(:,:)
    logical :: calculate_permutation_

    calculate_permutation_ = default_calculate_permutation
    if (present(calculate_permutation)) calculate_permutation_ = calculate_permutation
    
    A = A_input 
    k = min(size(A, dim=1), size(A, dim=2))
    call dgetrf2_wrapper(A, P(1 : k))

    call extract_LU_real_dp(A, L, U)

    if (calculate_permutation_) then
      row_pivot_map = P(1 : k)
      call pivot_to_permutation(row_pivot_map, P)
    end if
  end subroutine lu_row_pivot_real_dp

  !> See **[[lu_full_pivot_real_dp]]**. Use this interface if \( \mathbf{A} \) can be overwritten.
  !>
  !> This is a wrapper for the LAPACK routine dgetrf2.
  subroutine dgetrf2_wrapper(A, row_pivot_map)
    !> On entry, the matrix \( \mathbf{A} \), on exit 
    !> the factors \( \mathbf{L} \) and \( \mathbf{U} \), where the 
    !> unit diagonals of \( \mathbf{L} \) are not stored.
    real(dp), intent(inout), contiguous :: A(:, :)
    !> Pivot indices; for \( 1 \le \) `i` \(  \le N \), row  `A(i, :)`  has
    !> been interchanged with row `A(row_pivot_map(i), :)`.
    integer, intent(out), contiguous :: row_pivot_map(:)

    integer :: info
    
    call assert(size(row_pivot_map) == min(size(A, dim=1), size(A, dim=2)), &
               'The size of row_pivot_map must be equal to the minimum of rows and columns of A.')

    call dgetrf2(size(A, dim=1), size(A, dim=2), A, size(A, dim=1), row_pivot_map, info)

    call terminate_if_false(info >= 0, 'dgetrf2 failed.')
  end subroutine dgetrf2_wrapper


  !> Calculate the LU factorization with partial row pivoting for a complex matrix
  !> \( \mathbf{A} \in \mathbb{C}^{M \times N} \) such that
  !> \[
  !>   \mathbf{A} = \mathbf{P} \cdot \mathbf{L} \cdot \mathbf{U},
  !> \]
  !> where \( \mathbf{P} \) is the permutation matrix, \( \mathbf{L} \) is a lower triangular 
  !> (trapezoidal) matrix with unit diagonal elements and \( \mathbf{U} \) is upper triangular (trapezoidal).
  !> Use this interface, if  \( \mathbf{A} \) should be unchanged.
  !>
  !> This is a wrapper for the LAPACK routine zgetrf2.
  subroutine lu_row_pivot_complex_dp(A_input, L, U, P, calculate_permutation)
    !> Input matrix \( \mathbf{A} \in \mathbb{C}^{m \times n} \)
    complex(dp), intent(in), contiguous :: A_input(:, :)
    !> Lower triangular matrix  \( \mathbf{L} \) (or lower trapezoidal if \(m > n\)).
    complex(dp), intent(out), contiguous :: L(:, :)
    !> Upper triangular matrix \( \mathbf{U} \) (or upper trapezoidal if \(m < n\)).
    complex(dp), intent(out), contiguous :: U(:, :)
    !> If perm_or_piv is 'permutation', the permutation vector of the first min'(m, n)' rows is returned,
    !> else the  pivot vector of the first min'(m, n)' rows as obtained by **[[zgetrf2]]** is returned.
    integer, intent(out), contiguous :: P(:)
    !> Return either permutation vectors or pivot vectors.
    !> If `.true.`, the permutation vectors are returned,
    !> else the pivot vectors are returned.
    !> Default is `.false.`.
    logical, intent(in), optional :: calculate_permutation 
    
    integer :: k
    integer, allocatable :: row_pivot_map(:)
    complex(dp), allocatable :: A(:,:)
    logical :: calculate_permutation_

    calculate_permutation_ = default_calculate_permutation
    if (present(calculate_permutation)) calculate_permutation_ = calculate_permutation

    A = A_input
    k = min(size(A, dim=1), size(A, dim=2))
    call zgetrf2_wrapper(A, P(1 : k))

    call extract_LU_complex_dp(A, L, U)

    if (calculate_permutation_) then
      row_pivot_map = P(1 : k)
      call pivot_to_permutation(row_pivot_map, P)
    end if
  end subroutine lu_row_pivot_complex_dp


  !> See **[[lu_full_pivot_complex_dp]]**. Use this interface if \( \mathbf{A} \) can be overwritten.
  !>
  !> This is a wrapper for the LAPACK routine zgetrf2.
  subroutine zgetrf2_wrapper(A, row_pivot_map)
    !> On entry, the matrix \( \mathbf{A} \), on exit 
    !> the factors \( \mathbf{L} \) and \( \mathbf{U} \), where the 
    !> unit diagonals of \( \mathbf{L} \) are not stored.
    complex(dp), intent(inout), contiguous :: A(:, :)
    !> Pivot indices; for \( 1 \le \) `i` \(  \le N \), row  `A(i, :)`  has
    !> been interchanged with row `A(row_pivot_map(i), :)`.
    integer, intent(out), contiguous :: row_pivot_map(:)

    integer :: info
    
    call assert(size(row_pivot_map) == min(size(A, dim=1), size(A, dim=2)), &
               'The size of row_pivot_map must be equal to the minimum of rows and columns of A.')

    call zgetrf2(size(A, dim=1), size(A, dim=2), A, size(A, dim=1), row_pivot_map, info)

    call terminate_if_false(info >= 0, 'dgetrf2 failed.')
  end subroutine zgetrf2_wrapper

  !> For a LU factorization, calculated by **[[lu_row_pivot]],
  !> \[
  !>    \mathbf{A} = \mathbf{L} \cdot \mathbf{U},
  !> \]
  !> invert \( \mathbf{U} \) and solve the system of linear equations
  !> \[
  !>    \mathbf{A}^{-1} \cdot \mathbf{L} = \mathbf{U}^{-1}.
  !> \]
  !> This routine does not assert if the matrix \( \mathbf{A} \) is invertible.
  !>
  !> This is a wrapper for the LAPACK routine dgetri.
  subroutine dgetri_wrapper(LU, row_pivot_map)
    !> On entry, the factors \(\mathbf{L}\) and \(\mathbf{U}\) of \(\mathbf{A}\) 
    !> calculated by **[[lu_row_pivot]], on exit, the inverse of \(\mathbf{A}\).
    real(dp), intent(inout), contiguous :: LU(:, :)
    !> Pivot indices from **[[lu_row_pivot]]; 
    !> for \( 1 \le i \le N \), row \(i\) has been interchanged with row **row_pivot_map**(\(i\)).
    integer, intent(in), contiguous :: row_pivot_map(:)

    integer :: lwork, info
    real(dp), allocatable :: work(:)

    call assert(is_square(LU), 'LU is not a square matrix.')
    call assert(size(row_pivot_map) == size(LU, dim=1), 'row_pivot_map and LU have not the same size.')

    lwork = -1 
    allocate(work(1))
    call dgetri(size(LU, dim=1), LU, size(LU, dim=1), row_pivot_map, work, lwork, info)
    call terminate_if_false(info == 0, 'Workspace query for dgetri failed.')

    lwork = work(1)
    deallocate(work); allocate(work(lwork))
    call dgetri(size(LU, dim=1), LU, size(LU, dim=1), row_pivot_map, work, lwork, info)
    call terminate_if_false(info == 0, 'dgetri failed.')
  end subroutine dgetri_wrapper

  !> For a LU factorization, calculated by **[[lu_row_pivot]],
  !> \[
  !>    \mathbf{A} = \mathbf{L} \cdot \mathbf{U},
  !> \]
  !> invert \( \mathbf{U} \) and solve the system of linear equations
  !> \[
  !>    \mathbf{A}^{-1} \cdot \mathbf{L} = \mathbf{U}^{-1}.
  !> \]
  !> This routine does not assert if the matrix \( \mathbf{A} \) is invertible.
  !>
  !> This is a wrapper for the LAPACK routine zgetri.
  subroutine zgetri_wrapper(LU, row_pivot_map)
    !> On entry, the factors \(\mathbf{L}\) and \(\mathbf{U}\) of \(\mathbf{A}\) 
    !> calculated by **[[lu_row_pivot]], on exit, the inverse of \(\mathbf{A}\).
    complex(dp), intent(inout), contiguous :: LU(:, :)
    !> Pivot indices from **[[lu_row_pivot]]; 
    !> for \( 1 \le i \le N \), row \(i\) has been interchanged with row **row_pivot_map**(\(i\)).
    integer, intent(in), contiguous :: row_pivot_map(:)

    integer :: lwork, info
    complex(dp), allocatable :: work(:)

    call assert(is_square(LU), 'LU is not a square matrix.')
    call assert(size(row_pivot_map) == size(LU, dim=1), 'row_pivot_map and LU have not the same size.')

    lwork = -1 
    allocate(work(1))
    call zgetri(size(LU, dim=1), LU, size(LU, dim=1), row_pivot_map, work, lwork, info)
    call terminate_if_false(info == 0, 'Workspace query for dgetri failed.')

    lwork = work(1)
    deallocate(work); allocate(work(lwork))
    call zgetri(size(LU, dim=1), LU, size(LU, dim=1), row_pivot_map, work, lwork, info)
    call terminate_if_false(info == 0, 'dgetri failed.')
  end subroutine zgetri_wrapper

! utils

  !> Copy from a matrix \( \mathbf{A} \) the lower triangular part \( \mathbf{L} \) and the upper triangular part \( \mathbf{U} \) into 
  !> two arrays \( \mathbf{L} \) \( \mathbf{U} \), where the diagonal of \( \mathbf{L} \) is known to be 1.0.
  subroutine extract_LU_real_dp(A, L, U) 
    !> Input matrix \( \mathbf{A} \); contains 
    !> the factors \( \mathbf{L} \) and \( \mathbf{U} \), where the 
    !> unit diagonals of \( \mathbf{L} \) are not stored.
    real(dp), intent(in), contiguous :: A(:, :)
    !> Lower triangular matrix \( \mathbf{L} \)
    real(dp), intent(out), contiguous :: L(:, :)
    !> Upper triangular matrix \( \mathbf{U} \)
    real(dp), intent(out), contiguous :: U(:, :)

    integer :: i, m, n, k

    m = size(A, dim=1)
    n = size(A, dim=2)
    k = min(m, n)

! This preprocessor usage is deliberate, to prevent if statements being evaluated in production code.
#ifdef USE_ASSERT
    if (m > n) then
      call assert(all(shape(L) == shape(A)), 'If m > n, L must be of the same shape as A.')
      call assert(all(shape(U) == [k, k]), 'If m > n, U must be a suqare matrix of size min(m, n).')
    else 
      call assert(all(shape(L) == [k, k]), 'If m > n, L must be a suqare matrix of size min(m, n).')
      call assert(all(shape(U) == shape(A)), 'If m > n, U must be of the same shape as A.')
    end if 
#endif

    ! * Copy from A(1 : k, 1 : k-1) 
    !   - the lower triangular part to L(1 : k, 1 : k-1)
    !   - the upper triangular part to U(1 : k, 1 : k-1)
    ! * Set the diagonal of L to one
    ! * set the rest of L and U to zero
    do i=1, k-1
      L(1 : i-1, i) = 0.0_dp
      L(i, i) = 1.0_dp
      L(i+1 :, i) =  A(i+1 :, i)

      U(1 : i, i) = A(1 : i, i)
      U(i+1 :, i) = 0.0_dp
    end do

    ! Set the k'th column of L 
    L(:, k) = 0.0_dp
    L(k, k) = 1.0_dp
    ! Copy the k'th column of A to the k'th column of U
    U(: , k) = A(1 : k, k)

    ! If A has more rows than columns (m > n), copy A(k+1 : m, :)
    ! to L(k+1 : m, :)
    if (m > n) then
      L(k+1 : m, :) = A(k + 1 : m, :)
    
    ! If A has more columns than rows (m < n), copy A(:, k+1 : n)
    ! to U(:, k+1 : n)
    else if (m < n) then
      U(:, k+1 : n) = A(:, k+1 : n)
    end if 
  end subroutine extract_LU_real_dp


  !> Copy from a matrix \( \mathbf{A} \) the lower triangular part \( \mathbf{L} \) and the upper triangular part \( \mathbf{U} \) into 
  !> two arrays \( \mathbf{L} \) \( \mathbf{U} \), where the diagonal of \( \mathbf{L} \) is known to be (1.0, 0.0).
  subroutine extract_LU_complex_dp(A, L, U) 
    !> Input matrix \( \mathbf{A} \); contains 
    !> the factors \( \mathbf{L} \) and \( \mathbf{U} \), where the 
    !> unit diagonals of \( \mathbf{L} \) are not stored.
    complex(dp), intent(in), contiguous :: A(:, :)
    !> Lower triangular matrix \( \mathbf{L} \)
    complex(dp), intent(out), contiguous :: L(:, :)
    !> Upper triangular matrix \( \mathbf{U} \)
    complex(dp), intent(out), contiguous :: U(:, :)

    integer :: m, n, k

    m = size(A, dim=1)
    n = size(A, dim=2)
    k = min(m, n)

! This preprocessor usage is deliberate, to prevent if statements being evaluated in production code.
#ifdef USE_ASSERT
    if (m > n) then
      call assert(all(shape(L) == shape(A)), 'If m > n, L must be of the same shape as A.')
      call assert(all(shape(U) == [k, k]), 'If m > n, U must be a suqare matrix of size min(m, n).')
    else 
      call assert(all(shape(L) == [k, k]), 'If m > n, L must be a suqare matrix of size min(m, n).')
      call assert(all(shape(U) == shape(A)), 'If m > n, U must be of the same shape as A.')
    end if 
#endif

    ! * Copy from A(1 : k, 1 : k-1) 
    !   - the lower triangular part to L(1 : k, 1 : k-1)
    !   - the upper triangular part to U(1 : k, 1 : k-1)
    ! * Set the diagonal of L to one
    ! * set the rest of L and U to zero
    do i=1, k-1
      L(1 : i-1, i) = zzero
      L(i, i) = zone
      L(i+1 :, i) =  A(i+1 : m, i)

      U(1 : i, i) = A(1 : i, i)
      U(i+1 :, i) = zzero
    end do

    ! Set the k'th column of L 
    L(:, k) = zzero
    L(k, k) = zone
    ! Copy the k'th column of A to the k'th column of U
    U(:, k) = A(1 : k, k)
    
    ! If A has more rows than columns (m > n), copy A(k+1 : m, :)
    ! to L(k+1 : m, :)
    if (m > n) then
      L(k+1 : m, :) = A(k + 1 : m, :)
    
    ! If A has more columns than rows (m < n), copy A(:, k+1 : n)
    ! to U(:, k+1 : n)
    else if (m < n) then
      U(:, k+1 : n) = A(:, k+1 : n)
    end if 
  end subroutine extract_LU_complex_dp


  !> Calculate the permutation index map from the pivot index map. The permutation
  !> index map maps the index of the elements of an array to the index of the permuted array 
  !> such that `arr_permuted(i) = arr(permutation_map(i))`. The pivot index map maps the indices 
  !> to those which where interchanged, _e.g._, element `i` has been interchanged with element
  !> `pivot_map(i)`. To get the permutation index map from the pivot index map, one has to start
  !> with the identity permutation map (`[1, 2, 3, 4, ...]`) and redo the pivots for it, _e.g._ 
  !> \[ 
  !>    \begin{split}
  !>      \text{pivot_map}       &= (2, 2, 4, 1) \\\
  !>      \text{permutation_map} &= (1, 2, 3, 4)
  !>    \end{split}
  !> \]
  !> Step 1:
  !> Swap element 1 with element 2:
  !> \[ 
  !>    \begin{split}
  !>      \text{pivot_map}       &= (2, 2, 4, 1) \\\
  !>      \text{permutation_map} &= (2, 1, 3, 4)
  !>    \end{split}
  !> \]
  !> 
  !> Step 2:
  !> Swap element 2 with element 2:
  !> \[ 
  !>    \begin{split}
  !>      \text{pivot_map}       &= (2, 2, 4, 1) \\\
  !>      \text{permutation_map} &= (2, 1, 3, 4)
  !>    \end{split}
  !> \]
  !> 
  !> Step 3:
  !> Swap element 3 with element 4:
  !> \[ 
  !>    \begin{split}
  !>      \text{pivot_map}       &= (2, 2, 4, 1) \\\
  !>      \text{permutation_map} &= (2, 1, 4, 3)
  !>    \end{split}
  !> \]
  !> 
  !> Step 4:
  !> Swap element 4 with element 1:
  !> \[ 
  !>    \begin{split}
  !>      \text{pivot_map}       &= (2, 2, 4, 1) \\\
  !>      \text{permutation_map} &= (4, 2, 3, 1)
  !>    \end{split}
  !> \]
  !>
  !> Step \( i \):
  !> Swap element \( i \) of permutation_map with element pivot map \((i)\) of permutation_map:
  !> \[
  !>      \text{permutation map}(i) \leftrightarrow
  !>      \text{permutation map}( \text{pivot map} (i) )
  !> \]
  !>
  !> Where the last result of \(\text{permutation_map}\) is the permutation index map, according to the 
  !> pivot index map \( \text{pivot_map} \).
  subroutine pivot_to_permutation(pivot_map, permutation_map)
    !> Input array; contains pivot index map; for \( 1 \le i \le N \), 
    !> element \(i\) has been interchanged with element 'in(i)'.
    integer, intent(in), contiguous :: pivot_map(:)
    !> Output array; contains permutation indices such that
    !> the element at former position \(i\) is now at position 'out(i)'
    integer, intent(out), contiguous :: permutation_map(:)

    integer :: temp_pivot_index, temp_index

    call assert(size(pivot_map) <= size(permutation_map), &
            'The sizes of pivot_map must be less than or equal than the size of permutation_map.')
    
    permutation_map = mesh_1d(1, size(permutation_map))
    
    do i=1, size(pivot_map)
      temp_index = permutation_map(i)
      temp_pivot_index = permutation_map(pivot_map(i))

      permutation_map(i) = temp_pivot_index
      permutation_map(pivot_map(i)) = temp_index
    end do
  end subroutine pivot_to_permutation

end module lu_factorization