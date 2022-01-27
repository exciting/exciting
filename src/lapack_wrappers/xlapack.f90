!> Expose lapack wrapper
module xlapack
  ! multiplication
  use vector_multiplication, only: dot_multiply, outer_product
  use general_matrix_multiplication, only: matrix_multiply
  use hermitian_matrix_multiplication, only: hermitian_matrix_multiply
  ! decomposition
  use singular_value_decomposition, only: xgesdd, svd_divide_conquer
  ! utils
  use matrix_rank, only: matrix_rank_SVD
end module xlapack