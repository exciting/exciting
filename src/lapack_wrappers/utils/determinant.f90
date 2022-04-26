!> Calculate the determinant with different methods.
module determinant
  use asserts, only: assert
  use precision, only: dp
  use constants, only: zone
  use math_utils, only: is_square
  use lu_factorization, only: xgetc2
  use matrix_rank, only: matrix_rank_SVD

  private
  public :: determinant_LU

  !> Calculate the determinant of a matrix \( \mathbf{A} \in \mathbb{K}^{N \times N} \) by calculating the LU
  !> factorization with full pivoting. The determinant is than given by
  !> \[
  !>   \det(\mathbf{A}) = \det(P \cdot L \cdot U \cdot Q) = \det(P) \det(L) \det(U) \det(Q)
  !>           = (-1)^{S+T} \left( \prod_{i = 1}^{M} u_{ii} \right) \left( \prod_{i = 1}^{M} u_{ii} \right)
  !>           = (-1)^{S+T} \left( \prod_{i = 1}^{M} u_{ii} \right),
  !> \]
  !> where \( S, T \) are the numbers of row and column exchanges due to \( P, Q \) respectively and \( \mathbb{K} \)
  !> is one of \( \mathbb{R}, \mathbb{C} \).
  interface determinant_LU
    module procedure determinant_real_dp, &
                     determinant_complex_dp
  end interface determinant_LU

  contains

  !---------------------------------------------------------------------------------------------------------------------
  ! determinant_LU

  !> Calculate the determinant of a matrix \( \mathbf{A} \in \mathbb{R}^{N \times N} \) by calculating the LU
  !> factorization with full pivoting. (see **[[dgetc2_wrapper]]**). The determinant is than given by
  !> \[\begin{split}
  !>     \det(\mathbf{A}) &= \det(\mathbf{P} \cdot \mathbf{L} \cdot \mathbf{U} \cdot \mathbf{Q}) \\\
  !>             &= \det(\mathbf{P}) \det(\mathbf{L}) \det(\mathbf{U}) \det(\mathbf{Q}) \\\
  !>             &= (-1)^{S+T} \left( \prod_{i = 1}^{M} u_{ii} \right) \left( \prod_{i = 1}^{M} u_{ii} \right) \\\
  !>             &= (-1)^{S+T} \left( \prod_{i = 1}^{M} u_{ii} \right),
  !>   \end{split}
  !> \]
  !> where \( S, T \) are the numbers of row and column exchanges due to \( \mathbf{P}, \mathbf{Q} \) respectively.
  function determinant_real_dp(A_input) result(determinant)
    !> Input matrix
    real(dp), intent(in), contiguous :: A_input(:, :)

    real(dp) :: determinant

    integer :: S, T, i, n
    integer, allocatable :: row_pivot_map(:), column_pivot_map(:), indices_S(:), indices_T(:)
    logical, allocatable :: diff_S(:), diff_T(:)
    real(dp), allocatable :: A(:, :)

    call assert(is_square(A_input), 'A_input is not a square matrix.')
    n = size(A_input, dim=1)

    if (matrix_rank_SVD(A_input) < n) then
      determinant = 0.0_dp
      return
    end if

    A = A_input
    allocate(row_pivot_map(n))
    allocate(column_pivot_map(n))

    call xgetc2(A, row_pivot_map, column_pivot_map)

    allocate(diff_S(n))
    allocate(diff_T(n))
    determinant = 1._dp

    do i = 1, n
      determinant = determinant * A( i, i)
      diff_S(i) = i == row_pivot_map(i)
      diff_T(i) = i == column_pivot_map(i)
    end do

    indices_S = pack([(i, i=1, n)], diff_S)
    indices_T = pack([(i, i=1, n)], diff_T)

    S = size(indices_S)
    T = size(indices_T)

    determinant =  (-1)**(S + T) * determinant
  end function determinant_real_dp


  !> Calculate the determinant of a matrix \( \mathbf{A} \in \mathbb{C}^{N \times N} \) by calculating the LU
  !> factorization with full pivoting (see **[[zgetc2_wrapper]]**). The determinant is than given by
  !> \[\begin{split}
  !>     \det(\mathbf{A}) &= \det(\mathbf{P} \cdot \mathbf{L} \cdot \mathbf{U} \cdot \mathbf{Q}) \\\
  !>             &= \det(\mathbf{P}) \det(\mathbf{L}) \det(\mathbf{U}) \det(\mathbf{Q}) \\\
  !>             &= (-1)^{S+T} \left( \prod_{i = 1}^{M} u_{ii} \right) \left( \prod_{i = 1}^{M} u_{ii} \right) \\\
  !>             &= (-1)^{S+T} \left( \prod_{i = 1}^{M} u_{ii} \right),
  !>   \end{split}
  !> \]
  !> where \( S, T \) are the numbers of row and column exchanges due to \( \mathbf{P}, \mathbf{Q} \) respectively.
  function determinant_complex_dp(A_input) result(determinant)
    !> Input matrix
    complex(dp), intent(in), contiguous  :: A_input(:, :)

    complex(dp) :: determinant

    integer :: S, T, i, n
    integer, allocatable :: row_pivot_map(:), column_pivot_map(:), indices_S(:), indices_T(:)
    logical, allocatable :: diff_S(:), diff_T(:)
    complex(dp), allocatable :: A(:, :)

    call assert(is_square(A_input), 'A_input is not a square matrix.')
    n = size(A_input, dim=1)

    if (matrix_rank_SVD(A_input) < n) then
      determinant = 0.0_dp
      return
    end if

    A = A_input
    allocate(row_pivot_map(n))
    allocate(column_pivot_map(n))

    call xgetc2(A, row_pivot_map, column_pivot_map)

    allocate(diff_S(n))
    allocate(diff_T(n))
    determinant = zone

    do i = 1, n
      determinant = determinant * A( i, i)
      diff_S(i) = i == row_pivot_map(i)
      diff_T(i) = i == column_pivot_map(i)
    end do

    indices_S = pack([(i, i=1, n)], diff_S)
    indices_T = pack([(i, i=1, n)], diff_T)

    S = size(indices_S)
    T = size(indices_T)

    determinant =  (-1)**(S + T) * determinant
  end function determinant_complex_dp

end module determinant