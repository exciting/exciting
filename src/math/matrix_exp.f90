!> Module computing the matrix exponential
module matrix_exp
#include "asserts.fpp"
  use constants, only: zone, zzero
  use math_utils, only: is_hermitian, is_positive_definite
  use general_matrix_multiplication, only: matrix_multiply
  use hermitian_matrix_multiplication, only: hermitian_matrix_multiply
  use linear_system_positive_definite, only: positive_definite_solve
  use precision, only: dp, i32
  use xlapack, only: solve_generalized_hermitian_eigenproblem

  implicit none

  private
  public :: exp_hermitian_matrix_times_vectors, &
            exp_general_matrix_times_vectors, &
            exphouston_hermitian_matrix_times_vectors

  !> Default tolerance
  real(dp), parameter :: tol_default = 1e-6_dp

contains
  !> Obtain the action of the exponential \( \exp(\alpha \hat{H}) \) operator on
  !> a set of vectors, yielding \( \exp(\alpha \hat{H})| v_{j} \rangle \), 
  !> for \(j=1,\ldots, N\). Here, \( \hat{H} \) is a hermitian operator, and \( \alpha \)
  !> is a complex pre-factor. The vectors \( | v_{j} \rangle \)
  !> are represented in terms of expansion coefficients \( C_{j\mu} \)
  !> with respect to a non-orthonormal basis:
  !> \[
  !>    | v_{j} \rangle = \sum_\mu
  !>    C_{j\mu} | \phi_{\mu} \rangle
  !> \]
  !> Given the overlap matrix \( S \), the exponential operation is expressed as: 
  !> \[
  !>    \exp [ \alpha \hat{H} ] | v_{j} \rangle =
  !>    \exp [ \alpha S^{-1}H ] \; C_{j} = \sum_{n=0}^{M} \frac{1}{n!}
  !>    (\alpha S^{-1}H)^n C_{j}, \quad j=1,\ldots, N
  !> \]
  !> Here, the matrix exponential is approximated using a Taylor expansion
  !> up to the order defined by \( M \) ([[order_taylor]])
  subroutine exp_hermitian_matrix_times_vectors( order_taylor, alpha, &
    & H, S, vectors, tol )
    !> The order \( M \) of the Taylor expansion
    integer(i32), intent(in)       :: order_taylor
    !> Complex prefactor \( \alpha \)
    complex(dp), intent(in)        :: alpha
    !> Hermitian matrix \( H \)
    complex(dp),intent(in)         :: H(:, :)
    !> Overlap matrix \( S \): must be positive definite
    complex(dp),intent(in)         :: S(:, :)
    !> On entry: the expansion coefficients \( C_{j} \).
    !> On exit: \( \exp [ \alpha S^{-1}H ] \; C_{j}\).
    complex(dp),intent(inout)      :: vectors(:, :)
    !> Tolerance to check if matrices are hermitian and positive definite
    real(dp), intent(in), optional :: tol

    integer(i32)                  :: it, info, dim_H, n_vectors
    complex(dp), allocatable      :: x(:, :), y(:, :)
    real(dp)                      :: tolerance

    n_vectors = size( vectors, 2 )
    dim_H = size( H, 1 )
    allocate( x, source = vectors )
    allocate( y(dim_H, n_vectors) )

    ! Optional arguments
    tolerance = tol_default
    if( present(tol) ) tolerance = tol

    ! Sanity checks
    CALL_ASSERT( dim_H == size( vectors, 1 ), 'H and vectors have incompatible sizes.' )
    CALL_ASSERT( is_positive_definite( S, tolerance ), 'S is not positive definite' )
    CALL_ASSERT( size( S, 1 ) == size( vectors, 1 ), 'S and vectors have incompatible sizes.' )

    ! Taylor expansion
    do it = 1, order_taylor
      ! Matrix multiplication: y = H*x
      call hermitian_matrix_multiply( H, x, y, tol=tolerance )
      ! Obtain (S^(-1))*y for positive definite S (y will store the solution)
      call positive_definite_solve( S, y )
      x = ( alpha/it )*y
      vectors = vectors + x
    end do

  end subroutine 

  !> Same as [[exp_hermitian_matrix_times_vectors]], but for the case of 
  !> general matrix \(H\)
  subroutine exp_general_matrix_times_vectors( order_taylor, alpha, &
    & H, S, vectors, tol )
    !> The order of the Taylor expansion
    integer(i32), intent(in)      :: order_taylor
    !> Complex prefactor \( \alpha \)
    complex(dp), intent(in)       :: alpha
    !> General matrix \( H \)
    complex(dp),intent(in)        :: H(:, :)
    !> Overlap matrix \( S \): must be positive definite
    complex(dp),intent(in)        :: S(:, :)
    !> On exit: \( \exp [ \alpha S^{-1}H ] C\)
    complex(dp),intent(inout)     :: vectors(:, :)
    !> Tolerance to check if matrices are hermitian and positive definite
    real(dp), intent(in), optional:: tol

    integer(i32)                  :: it, info
    integer(i32)                  :: dim_H, n_vectors
    complex(dp), allocatable      :: x(:, :), y(:, :)
    real(dp)                      :: tolerance

    ! Allocate arrays
    n_vectors = size( vectors, 2 )
    dim_H = size( H, 1 )
    allocate( x, source = vectors )
    allocate( y(dim_H, n_vectors) )

    ! Optional arguments
    tolerance = tol_default
    if( present(tol) ) tolerance = tol

    ! Sanity checks
    CALL_ASSERT( size( H, 1 ) == size( vectors, 1 ), 'H and vectors have incompatible sizes.' )
    CALL_ASSERT( is_positive_definite( S, tolerance ), 'S is not positive definite' )
    CALL_ASSERT( size( S, 1 ) == size( vectors, 1 ), 'S and vectors have incompatible sizes.' )

    ! Taylor expansion
    do it = 1, order_taylor
      ! Matrix multiplication: y = H*x
      call matrix_multiply( H, x, y )
      ! Obtain (S^(-1))*y for positive definite S (y will store the solution)
      call positive_definite_solve( S, y )
      x = ( alpha/it )*y
      vectors = vectors + x
    end do

  end subroutine

  !> Similar to [[exp_hermitian_matrix_times_vectors]], but without employing a Taylor expansion.
  !> Instead, here the Houston expansion (see this
  !> [paper](https://doi.org/10.1103/PhysRevB.89.224305)) is used, which allows for an
  !> exact evaluation of the exponential operator. This is achieved by introducing an 
  !> auxiliary basis composed of the eigenvectors of \( \hat{H} \). 
  !> Specifically, the eigenvalue problem
  !> \( \hat{H}| u^0_{i}\rangle = \varepsilon_{i}| u^0_{i}\rangle\) is solved, where
  !> \( i \) runs from 1 to \( N \leq \) \( \dim( \hat{H} ) \). Then:
  !> \[ 
  !>     \exp(\alpha \hat{H}) = \sum_{i = 1}^{N} 
  !>     | u^0_i\rangle \exp(\alpha \varepsilon_i) \langle u^0_i |.
  !> \]
  !> In practice, the generalized eigenvalue problem is solved:
  !>  \[
  !>       H C^0_{i} = \varepsilon_{i} S C^0_{i},
  !>  \]
  !> where \( C^0_{i} \) represents the expansion coefficients of 
  !> \( | u^0_{i}\rangle \) in terms of a non-orthonormal basis
  !> \( | \phi_{\mu} \rangle \) with overlap matrix \(S\).
  !> \[  
  !>      | u^0_{j} \rangle = \sum_\mu  C^0_{i\mu} | \phi_{\mu} \rangle. 
  !> \]
  !> The goal is to express an arbitrary state \( | v_{j}\rangle \) in terms of
  !> the eigenstates \( | u^0_{i}\rangle \)
  !>  \[
  !>      | v_{j}\rangle = \sum_i p_{ij}
  !>      | u^0_{i}\rangle,
  !>  \]
  !> where the projection coefficients are given by \( p_{ij} =
  !> \langle u^0_{i} | v_{j}\rangle \):
  !> \[
  !>     p_{ij} = (C^0_{i})^\dagger S C_{j},
  !> \]
  !> Since the action of the exponential operator on an eigenstate is given by:  
  !>  \[
  !>	   \exp(\alpha \hat{H}) |u^0_{i}\rangle =
  !>     \exp [ \alpha S^{-1}H ] C^0_{i} =
  !>     \mathrm{e}^{\alpha\varepsilon_{i}} C^0_{i},
  !>  \]
  !> the desired operation follows as
  !>  \[
  !>        \exp(\alpha \hat{H})	| v_{j}\rangle
  !>          =   \sum_i p_{ij} \exp(\alpha \hat{H})
  !>        | u^0_{i}\rangle =
  !>       \sum_i  \mathrm{e}^{\alpha\varepsilon_{i}}
  !>      C^0_{i} p_{ij} =
  !>      \sum_i \tilde{C}^0_{i} p_{ij},
  !> \]
  !> where \( \tilde{C}^0_{i} = \mathrm{e}^{\alpha\varepsilon_{i}}  C^0_{i} \).
  subroutine exphouston_hermitian_matrix_times_vectors( alpha, H, S, vectors, tol, n_eigs )
    !> Complex pre-factor \( \alpha \)
    complex(dp), intent(in) :: alpha
    !> Hermitian matrix \( H \)
    complex(dp),intent(in) :: H(:, :)
    !> Overlap matrix \( S \): must be positive definite
    complex(dp),intent(in) :: S(:, :)
    !> Refer to [[exp_hermitian_matrix_times_vectors]]
    complex(dp),intent(inout) :: vectors(:, :)
    !> Tolerance for checking if the matrices are hermitian / positive definite
    real(dp), optional, intent(in) :: tol
    !> Number \( N \) of the lowest-lying eigenstates used for the operator expansion
    integer(i32), optional, intent(in) :: n_eigs

    integer(i32) :: i, n_vectors, matrix_dim, n_expansion
    real(dp) :: tolerance
    real(dp), allocatable :: eigenvalues(:)
    complex(dp), allocatable :: eigenvectors(:, :), proj(:, :), aux(:, :), aux_exp(:, :)
    complex(dp), allocatable :: S_copy(:, :), H_copy(:, :)

    tolerance = tol_default
    if( present(tol) ) tolerance = tol
    matrix_dim = size( H, 1 )
    n_expansion = matrix_dim
    if ( present( n_eigs ) ) n_expansion = n_eigs
    n_vectors = size( vectors, 2 )

    ! Sanity checks
    CALL_ASSERT( n_expansion <= matrix_dim, 'more eigenvectors than matrix_dim are requested' )
    CALL_ASSERT( size( H, 1 ) == size( vectors, 1 ), 'H and vectors have incompatible sizes.' )
    CALL_ASSERT( is_positive_definite( S, tolerance ), 'S is not positive definite' )
    CALL_ASSERT( size( S, 1 ) == size( vectors, 1 ), 'S and vectors have incompatible sizes.' )

    allocate( eigenvectors(matrix_dim, n_expansion), aux(matrix_dim, n_vectors) )
    allocate( proj(n_expansion, n_vectors), aux_exp(matrix_dim, n_expansion), eigenvalues(n_expansion) )
    allocate( S_copy, source = S ) 
    allocate( H_copy, source = H )

    call solve_generalized_hermitian_eigenproblem( H_copy, S_copy, tol, eigenvalues, eigenvectors )

    ! Project vectors onto the eigenvectors
    call hermitian_matrix_multiply( S, vectors, aux, tol=tolerance )
    call matrix_multiply( eigenvectors, aux, proj, 'C')

    ! Scale each eigenvector by the exponential of alpha*eigenvalues
    do concurrent (i = 1: n_expansion)
      aux_exp(:, i) = exp( alpha * eigenvalues(i) ) * eigenvectors(:, i)
    end do

    call matrix_multiply( aux_exp, proj, vectors )

  end subroutine 


end module matrix_exp
