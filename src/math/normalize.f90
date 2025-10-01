
module normalize

  use asserts, only: assert
  use precision, only: dp, i32
  use math_utils, only: is_positive_definite
  use xlapack, only: dot_multiply, hermitian_matrix_multiply
  implicit none

  private 
  public :: normalize_vectors, &
            norm_squared_with_positive_matrix

  interface normalize_vectors
    module procedure :: normalize_vectors_with_matrix_complex_dp
  end interface normalize_vectors

  !> Obtain the norm squared of each vector in a set. The norm is calculated using 
  !> a positive definite matrix `S`.
  interface norm_squared_with_positive_matrix
    module procedure :: norm_squared_with_positive_matrix_complex_dp
  end interface norm_squared_with_positive_matrix


contains

  !> Normalize a set of vectors \( x_j \) dividing each vector by its norm as defined 
  !> in [[norm_squared_with_positive_matrix]]
  subroutine normalize_vectors_with_matrix_complex_dp( S, vectors )
    !> Overlap matrix: must be positive definite.
    complex(dp), contiguous, intent(in)        :: S(:, :)
    !> Vectors to be normalized stored as columns.
    complex(dp), contiguous, intent(inout)     :: vectors(:, :)

    integer(i32) :: j
    real(dp), allocatable :: norms_squared(:)

    associate( n_vectors => size(vectors, 2) )
      allocate( norms_squared(n_vectors) )
      call norm_squared_with_positive_matrix( vectors, S, norms_squared )
      do j = 1, n_vectors
        vectors(:, j) = vectors(:, j)/sqrt( norms_squared(j) ) 
      end do
    end associate
  end subroutine

  !> Obtain the norm squared of each vector in a set as
  !> \[
  !>  x^\dagger S x,
  !> \]
  !> where \( S \) is a positive definite matrix.  
  !> Note that this is the dot product between \( x \) and 
  !> \( S x \).
  subroutine norm_squared_with_positive_matrix_complex_dp( vectors, matrix, norms_squared )
    !> Vectors stored as columns to obtain their norm squared
    complex(dp), contiguous, intent(in) :: vectors(:, :)
    !> Positive definite matrix `S`
    complex(dp), contiguous, intent(in) :: matrix(:, :)
    !> The norm squared of each vector in `vectors`
    real(dp), contiguous, intent(out) :: norms_squared(:)

    complex(dp), allocatable :: A(:, :)
    integer(i32) :: j

    call assert( is_positive_definite(matrix), 'matrix is not positive definite.' )
    associate( m => size(matrix, 1), n_vectors => size(vectors, 2) )
      call assert( size(norms_squared) == n_vectors, 'norms must have size = n_vectors' )
      allocate( A(m, n_vectors) )
      call hermitian_matrix_multiply( matrix, vectors, A )
      do j = 1, n_vectors
        norms_squared(j) = real( dot_multiply( vectors(:, j), A(:, j), conjg_a = .true. ), dp )
      end do
    end associate
  end subroutine
end module