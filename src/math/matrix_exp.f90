!> Module for advanced matrix operations
module matrix_exp
  use asserts, only: assert
  use constants, only: zone, zzero
  use math_utils, only: is_hermitian, is_positive_definite
  use general_matrix_multiplication, only: matrix_multiply
  use hermitian_matrix_multiplication, only: hermitian_matrix_multiply
  use precision, only: dp

  implicit none

  private
  public :: exp_hermitianoperator_times_wavefunctions, &
            exphouston_hermitianoperator_times_wavefunctions

contains
  !> This subroutine obtains the exponential \( \exp(\alpha \hat{H}) \) applied
  !> to a set of vectors: \( \exp(\alpha \hat{H})| \Psi_{j\mathbf{k}} \rangle \).
  !> The operator \( \hat{H} \) **must** be hermitian. \( \alpha \)
  !> is a complex prefactor. The vectors \( | \Psi_{j\mathbf{k}} \rangle \)
  !> are described through the expansion coefficients \( C_{j\mathbf{k}\mu} \)
  !> in terms of the (L)APW+lo basis.
  !> \( | \phi_{\mathbf{k}\mu} \rangle \)
  !> \[
  !>    | \Psi_{j\mathbf{k}} \rangle = \sum_\mu
  !>    C_{j\mathbf{k}\mu} | \phi_{\mathbf{k}\mu} \rangle
  !> \]
  !> Since the basis is not orthonormal, we have
  !> \[
  !>    \exp [ \alpha \hat{H} ]
  !>    | \Psi_{j\mathbf{k}} \rangle =
  !>    \exp [ \alpha S_{\mathbf{k}}^{-1}H_{\mathbf{k}} ] \;
  !>    C_{j\mathbf{k}} = \sum_{n=0}^{M} \frac{1}{n!}
  !>    (\alpha S_{\mathbf{k}}^{-1}H_{\mathbf{k}})^n C_{j\mathbf{k}}
  !> \]
  !> The exponential here is approximated by a Taylor expansion
  !> up to the order defined by \( M \) (`order_taylor`)
  subroutine exp_hermitianoperator_times_wavefunctions( order_taylor, alpha, &
    & H, S, vectors )
    !> The order of the Taylor expansion
    integer, intent(in)           :: order_taylor
    !> Complex prefactor
    complex(dp), intent(in)       :: alpha
    !> Hermitian matrix \( H_{\mathbf{k}} \)
    complex(dp),intent(in)        :: H(:, :)
    !> Overlap matrix \( S_{\mathbf{k}} \)
    complex(dp),intent(in)        :: S(:, :)
    !> On entry: the expansion coefficients of
    !> \( | \Psi_{j\mathbf{k}} \rangle \) in terms of (L)APW+lo.
    !> On exit: \( \exp [ \alpha S_{\mathbf{k}}^{-1}H_{\mathbf{k}} ] \;
    !>    C_{j\mathbf{k}}\)
    complex(dp),intent(inout)     :: vectors(:, :)

    integer                       :: it, info
    integer                       :: dim, n_vectors
    complex(dp), allocatable      :: x(:, :), y(:, :), S_copy(:, :)



    ! Allocate arrays
    n_vectors = size( vectors, 2 )
    dim = size( H, 1 )
    allocate( x(dim, n_vectors) )
    allocate( y(dim, n_vectors) )
    allocate( S_copy(dim, dim) )
    ! Sanity checks
    ! Check if H is hermitian
    call assert( is_hermitian( H ), 'H is not hermitian' )
    ! Check if H and vectors have compatible size
    call assert( size( H, 1 ) == size( vectors, 1 ), 'H and vectors have incompatible sizes.' )
    ! Check if S is positive definite
    call assert( is_positive_definite( S ), 'S is not positive definite' )
    ! Check if S and vectors have compatible size
    call assert( size( S, 1 ) == size( vectors, 1 ), 'S and vectors have incompatible sizes.' )

    ! Initializations
    x = vectors
    S_copy = S

    ! Taylor expansion
    do it = 1, order_taylor
      ! Matrix multiplication: y = H*x
      call hermitian_matrix_multiply( H, x, y )
      ! Obtain (S^(-1))*y for positive definite S (y will store the solution)
      call ZPOSV( 'U', dim, n_vectors, S_copy, dim, y, dim, info )
      ! Restores S_copy to its original value, after being modified by ZPOSV
      S_copy = S
      x = ( alpha/it )*y
      vectors = vectors + x
    end do

  end subroutine exp_hermitianoperator_times_wavefunctions


  !> This subroutine obtains the exponential \( \exp(\alpha \hat{H}) \) applied
  !> to a set of vectors \( | \Psi_{j\mathbf{k}}\rangle \),
  !> as done in [[exp_hermitianoperator_times_wavefunctions]].
  !> Here the so-called Houston expansion (see this
  !> [paper](https://doi.org/10.1103/PhysRevB.89.224305)) is employed.
  !> We evaluate the exponential operator exactly rather than Taylor-expanding
  !> it. This can be done by taking into account an auxiliary basis formed by
  !> the eigenvectors of \( \hat{H} \). This means that we solve
  !> \( \hat{H}| \psi^0_{i\mathbf{k}}\rangle =
  !>       \varepsilon_{i\mathbf{k}}| \psi^0_{i\mathbf{k}}\rangle\), which in
  !> practice is carried out solving:
  !>  \[
  !>       H_{\mathbf{k}} C^0_{i\mathbf{k}} =
  !>       \varepsilon_{i\mathbf{k}} S_{\mathbf{k}} C^0_{i\mathbf{k}},
  !>  \]
  !> where \( C^0_{i\mathbf{k}} \) is an array to represent the expansions
  !> coefficients of \( | \psi^0_{i\mathbf{k}}\rangle \) in terms of the basis
  !> \( | \phi_{\mathbf{k}\mu} \rangle \)
  !>  \[  | \psi^0_{j\mathbf{k}} \rangle = \sum_\mu
  !>      C^0_{i\mathbf{k}\mu} | \phi_{\mathbf{k}\mu}. \rangle \]
  !> We now want to write \( | \Psi_{j\mathbf{k}}\rangle \) in terms of
  !> \( | \psi^0_{i\mathbf{k}}\rangle \)
  !>  \[
  !>      | \Psi_{j\mathbf{k}}\rangle = \sum_i p_{ij\mathbf{k}}
  !>      | \psi^0_{i\mathbf{k}}\rangle,
  !>  \]
  !> where the projection coefficients are given by \( p_{ij\mathbf{k}} =
  !> \langle \psi^0_{i\mathbf{k}} | \Psi_{j\mathbf{k}}\rangle \).
  !> If \( | \Psi_{j\mathbf{k}}\rangle \) are represented by their expansion
  !> coefficients \( C_{j\mathbf{k}} \), then we calculate
  !> \( p_{ij\mathbf{k}} \) through
  !> \[
  !>     p_{ij\mathbf{k}} = (C^0_{i\mathbf{k}})^\dagger S_{\mathbf{k}}
  !>                  C_{j\mathbf{k}}.
  !> \]
  !> Since
  !>  \[
  !>	   \exp(\alpha \hat{H}) |\psi^0_{i\mathbf{k}}\rangle =
  !>     \exp [ \alpha S_{\mathbf{k}}^{-1}H_{\mathbf{k}} ] C^0_{i\mathbf{k}} =
  !>     \mathrm{e}^{\alpha\varepsilon_{i\mathbf{k}}}
  !>     C^0_{i\mathbf{k}},
  !>  \]
  !> therefore, coming back to our original task, we have
  !>  \[
  !>        \exp(\alpha \hat{H})	| \Psi_{j\mathbf{k}}\rangle
  !>          =   \sum_i p_{ij\mathbf{k}} \exp(\alpha \hat{H})
  !>        | \psi^0_{i\mathbf{k}}\rangle =
  !>       \sum_i  \mathrm{e}^{\alpha\varepsilon_{i\mathbf{k}}}
  !>      C^0_{i\mathbf{k}} p_{ij\mathbf{k}} =
  !>      \sum_i \tilde{C}^0_{i\mathbf{k}} p_{ij\mathbf{k}},
  !> \]
  !> where \( \tilde{C}^0_{i\mathbf{k}} =
  !>          \mathrm{e}^{\alpha\varepsilon_{i\mathbf{k}}}
  !>           C^0_{i\mathbf{k}}\).
  subroutine exphouston_hermitianoperator_times_wavefunctions( alpha, &
      & H, S, vectors, tol )
    !> Complex prefactor
    complex(dp), intent(in)   :: alpha
    !> Hermitian matrix \( H_{\mathbf{k}} \)
    complex(dp),intent(in)    :: H(:, :)
    !> Overlap matrix
    complex(dp),intent(in)    :: S(:, :)
    !> Refer to [[exp_hermitianmatrix_times_vectors]]
    complex(dp),intent(inout) :: vectors(:, :)
    !> Tolerance for the eigenvalue problem
    real(dp), intent(in), optional :: tol

    integer                   :: i, lwork, info, n_eigvals_found
    integer                   :: n_vectors, dim
    integer, allocatable      :: ifail(:), iwork(:)
    real(dp)                  :: tolerance
    real(dp)                  :: vl, vu
    real(dp), allocatable     :: eigvals(:)
    complex(dp), allocatable  :: rwork(:), work(:)
    complex(dp), allocatable  :: eigvecs(:, :), proj(:, :), aux(:, :)
    complex(dp), allocatable  :: S_copy(:, :), H_copy(:, :)

    tolerance = 1.0e-8_dp
    if( present(tol) ) tolerance = tol
    dim = size( H, 1 )
    n_vectors = size( vectors, 2 )
    allocate( eigvecs(dim, n_vectors), aux(dim, n_vectors) )
    allocate( proj(n_vectors, n_vectors) )
    allocate( S_copy(dim, dim), H_copy(dim, dim) )
    allocate( ifail(dim), iwork(5*dim), eigvals(dim), rwork(7*dim) )

    dim = size( H, 1 )
    ! Sanity checks
    ! Check if H is hermitian
    call assert( is_hermitian( H ), 'H is not hermitian' )
    ! Check if H and vectors have compatible size
    call assert( size( H, 1 ) == size( vectors, 1 ), 'H and vectors have incompatible sizes.' )
    ! Check if S is positive definite
    call assert( is_positive_definite( S ), 'S is not positive definite' )
    ! Check if S and vectors have compatible size
    call assert( size( S, 1 ) == size( vectors, 1 ), 'S and vectors have incompatible sizes.' )

    S_copy = S
    H_copy = H
    vl = 0._dp
    vu = 0._dp

    ! Obtain the optimum lwork
    lwork = -1
    allocate( work(2) )
    call ZHEGVX( 1, 'V', 'I', 'U', dim, H_copy, dim, S_copy, dim, vl, &
      & vu, 1, n_vectors, tol, n_eigvals_found, eigvals, eigvecs, dim, &
      & work, lwork, rwork, iwork, ifail, info )
    lwork = int( work(1) )
    deallocate( work )
    allocate( work(lwork) )

    ! Solve the generalized eigenvalue/eigenvector problem: H*x = lambda*S*x
    ! We use H_copy and S_copy because ZHEGVX overwrites these matrices
    call ZHEGVX( 1, 'V', 'I', 'U', dim, H_copy, dim, S_copy, dim, vl, &
      & vu, 1, n_vectors, tol, n_eigvals_found, eigvals, eigvecs, dim, &
      & work, lwork, rwork, iwork, ifail, info )

    ! Check if there were problems with the diagonalization
    call assert( info==0, 'exphouston_hermitianmatrix_times_vectors: problems &
      & with ZHEGVX, info not zero' )

    ! Project vectors onto the eigenvectors
    call hermitian_matrix_multiply( S, vectors, aux )
    call matrix_multiply( eigvecs, aux, proj, 'C')

    ! Now, scale each eigenvector by the exponential of alpha*eigvals
    forall( i = 1:n_vectors ) aux(:, i) = zexp( alpha*eigvals(i) )*eigvecs(:, i)

    call matrix_multiply( aux, proj, vectors )

  end subroutine exphouston_hermitianoperator_times_wavefunctions


end module matrix_exp