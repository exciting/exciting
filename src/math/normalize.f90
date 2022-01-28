
module normalize

  use asserts, only: assert
  use precision, only: dp
  use math_utils, only: is_positive_definite
  use hermitian_matrix_multiplication, only: hermitian_matrix_multiply
  implicit none

  private 
  public :: normalize_vectors

  interface normalize_vectors
    module procedure :: normalize_vectors_with_matrix_complex_dp
  end interface normalize_vectors


contains

  !> Normalize a set of vectors. For each vector, the norm squared is
  !> \[
  !>  C_{j\mathbf{k}}^\dagger S_{\mathbf{k}} C_{j\mathbf{k}}
  !>  \]
  !> where \( C_{j\mathbf{k}} \) is the array of coefficients of the
  !> wavefunction in terms of the (L)APW+lo basis, and \( S_{\mathbf{k}} \)
  !> is the overlap matrix.  
  !> Note that, for each \(j\),
  !> \( C_{j\mathbf{k}}^\dagger S_{\mathbf{k}} C_{j\mathbf{k}} \) is the dot
  !> product between \( C_{j\mathbf{k}} \) and 
  !> \( S_{\mathbf{k}} C_{j\mathbf{k}} \).
  subroutine normalize_vectors_with_matrix_complex_dp( S, vectors )
    !> Overlap matrix: must be positive definite.
    complex(dp),intent(in)        :: S(:, :)
    !> Vectors to be normalized.
    !> First dimension: number of basis elements. Must be equal to each
    !> dimension of `S`.
    !> Second dimension: number of vectors
    complex(dp),intent(inout)     :: vectors(:, :)

    integer                       :: j, dim
    integer                       :: n_vectors
    real(dp)                      :: norm
    complex(dp), allocatable      :: A(:, :)
    complex(dp), external         :: zdotc

    call assert( is_positive_definite(S), 'S is not positive definite.' )

    n_vectors = size( vectors, 2 )
    dim = size( S, 1 )
    allocate( A(dim, n_vectors) )

    ! Matrix multiplication: A = S*vectors
    call hermitian_matrix_multiply( S, vectors, A )
    do j = 1, n_vectors
      ! Dot product between the j-th vector and the j-th column of S*vectors = A
      norm = dsqrt( dble( zdotc( dim, vectors(:, j), 1, A(:, j), 1 ) ) )
      vectors(:, j) = vectors(:, j)/norm
    end do
  end subroutine
end module