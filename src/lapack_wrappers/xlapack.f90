!> Expose lapack wrapper
module xlapack
  ! multiplication
  use vector_multiplication, only: dot_multiply
  use general_matrix_multiplication, only: matrix_multiply
  use hermitian_matrix_multiplication, only: hermitian_matrix_multiply
end module xlapack