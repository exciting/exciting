!> Expose lapack wrapper
module xlapack
  ! lapack workspaces
  use lapack_workspaces, only: lapack_workspace_complex_dp_t
  ! multiplication
  use add_scaled_array, only: scaled_add
  use vector_multiplication, only: dot_multiply, norm, outer_product
  use general_matrix_multiplication, only: matrix_multiply
  use hermitian_matrix_multiplication, only: hermitian_matrix_multiply
  ! decomposition
  use singular_value_decomposition, only: xgesdd, svd_divide_conquer
  use lu_factorization, only: lu_full_pivot, lu_row_pivot, xgetc2, xgetrf2, xgetri
  use qr_factorization, only: qr_column_pivot, xgeqp3, xorgqr, extract_R
  ! diagonalization
  use diagonalize_tridiagonal, only: diagonalize_symtridiag
  use generalized_hermitian_eigenproblem, only: solve_generalized_hermitian_eigenproblem
  ! linear solvers
  use linear_system_ill_defined_safe, only: ill_defined_safe_solve
  ! utils
  use matrix_rank, only: matrix_rank_SVD
  use determinant, only: determinant_LU
  use inverse, only: inverse_LU, invert_LU
end module xlapack
