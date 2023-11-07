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
  public :: exp_hermitian_matrix_times_vectors, &
            exp_general_matrix_times_vectors, &
            exphouston_hermitian_matrix_times_vectors

  !> Default tolerance
  real(dp), parameter :: tol_default = 1e-6_dp

contains
  !> This subroutine obtains the exponential \( \exp(\alpha \hat{H}) \) applied
  !> to a set of vectors: \( \exp(\alpha \hat{H})| v_{j} \rangle \).
  !> The operator \( \hat{H} \) **must** be hermitian. \( \alpha \)
  !> is a complex prefactor. The vectors \( | v_{j} \rangle \)
  !> are described through the expansion coefficients \( C_{j\mu} \)
  !> in terms of the (L)APW+lo basis.
  !> \[
  !>    | v_{j} \rangle = \sum_\mu
  !>    C_{j\mu} | \phi_{\mu} \rangle
  !> \]
  !> Since the basis is not orthonormal, we have
  !> \[
  !>    \exp [ \alpha \hat{H} ]
  !>    | v_{j} \rangle =
  !>    \exp [ \alpha S^{-1}H ] \;
  !>    C_{j} = \sum_{n=0}^{M} \frac{1}{n!}
  !>    (\alpha S^{-1}H)^n C_{j}
  !> \]
  !> The exponential here is approximated by a Taylor expansion
  !> up to the order defined by \( M \) (`order_taylor`)
  subroutine exp_hermitian_matrix_times_vectors( order_taylor, alpha, &
    & H, S, vectors, tol )
    !> The order of the Taylor expansion
    integer, intent(in)           :: order_taylor
    !> Complex prefactor
    complex(dp), intent(in)       :: alpha
    !> Hermitian matrix \( H \)
    complex(dp),intent(in)        :: H(:, :)
    !> Overlap matrix \( S \): must be positive definite
    complex(dp),intent(in)        :: S(:, :)
    !> On entry: the expansion coefficients of
    !> \( | v_{j} \rangle \) in terms of (L)APW+lo.
    !> On exit: \( \exp [ \alpha S_{\mathbf{k}}^{-1}H_{\mathbf{k}} ] \;
    !>    C_{j}\)
    complex(dp),intent(inout)     :: vectors(:, :)
    !> Tolerance to check if matrices are hermitian and positive definite
    real(dp), intent(in), optional:: tol
    integer                       :: it, info
    integer                       :: dim, n_vectors
    complex(dp), allocatable      :: x(:, :), y(:, :), S_copy(:, :)
    real(dp)                      :: tolerance



    ! Allocate arrays
    n_vectors = size( vectors, 2 )
    dim = size( H, 1 )
    allocate( x, source = vectors )
    allocate( y(dim, n_vectors) )
    allocate( S_copy, source=S )

    ! Optional arguments
    tolerance = tol_default
    if( present(tol) ) tolerance = tol

    ! Sanity checks
    ! Check if H is hermitian
    call assert( is_hermitian( H, tolerance ), 'H is not hermitian' )
    ! Check if H and vectors have compatible size
    call assert( size( H, 1 ) == size( vectors, 1 ), 'H and vectors have incompatible sizes.' )
    ! Check if S is positive definite
    call assert( is_positive_definite( S, tolerance ), 'S is not positive definite' )
    ! Check if S and vectors have compatible size
    call assert( size( S, 1 ) == size( vectors, 1 ), 'S and vectors have incompatible sizes.' )

    ! Taylor expansion
    do it = 1, order_taylor
      ! Matrix multiplication: y = H*x
      call hermitian_matrix_multiply( H, x, y, tol=tolerance )
      ! Obtain (S^(-1))*y for positive definite S (y will store the solution)
      call ZPOSV( 'U', dim, n_vectors, S_copy, dim, y, dim, info )
      ! Restores S_copy to its original value, after being modified by ZPOSV
      S_copy = S
      x = ( alpha/it )*y
      vectors = vectors + x
    end do

  end subroutine 

  !> Same as [[exp_hermitian_matrix_times_vectors]], but for the case of 
  !> general matrix \(H\)
  subroutine exp_general_matrix_times_vectors( order_taylor, alpha, &
    & H, S, vectors, tol )
    !> The order of the Taylor expansion
    integer, intent(in)           :: order_taylor
    !> Complex prefactor
    complex(dp), intent(in)       :: alpha
    !> General matrix \( H \)
    complex(dp),intent(in)        :: H(:, :)
    !> Overlap matrix \( S \): must be positive definite
    complex(dp),intent(in)        :: S(:, :)
    !> On exit: \( \exp [ \alpha S^{-1}H ] C\)
    complex(dp),intent(inout)     :: vectors(:, :)
    !> Tolerance to check if matrices are hermitian and positive definite
    real(dp), intent(in), optional:: tol
    integer                       :: it, info
    integer                       :: dim, n_vectors
    complex(dp), allocatable      :: x(:, :), y(:, :), S_copy(:, :)
    real(dp)                      :: tolerance

    ! Allocate arrays
    n_vectors = size( vectors, 2 )
    dim = size( H, 1 )
    allocate( x, source = vectors )
    allocate( y(dim, n_vectors) )
    allocate( S_copy, source=S )

    ! Optional arguments
    tolerance = tol_default
    if( present(tol) ) tolerance = tol

    ! Sanity checks
    ! Check if H and vectors have compatible size
    call assert( size( H, 1 ) == size( vectors, 1 ), 'H and vectors have incompatible sizes.' )
    ! Check if S is positive definite
    call assert( is_positive_definite( S, tolerance ), 'S is not positive definite' )
    ! Check if S and vectors have compatible size
    call assert( size( S, 1 ) == size( vectors, 1 ), 'S and vectors have incompatible sizes.' )

    ! Taylor expansion
    do it = 1, order_taylor
      ! Matrix multiplication: y = H*x
      call matrix_multiply( H, x, y )
      ! Obtain (S^(-1))*y for positive definite S (y will store the solution)
      call ZPOSV( 'U', dim, n_vectors, S_copy, dim, y, dim, info )
      ! Restores S_copy to its original value, after being modified by ZPOSV
      S_copy = S
      x = ( alpha/it )*y
      vectors = vectors + x
    end do

  end subroutine

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
  subroutine exphouston_hermitian_matrix_times_vectors( alpha, &
      & H, S, vectors, tol )
    !> Complex prefactor
    complex(dp), intent(in)   :: alpha
    !> Hermitian matrix \( H_{\mathbf{k}} \)
    complex(dp),intent(in)    :: H(:, :)
    !> Overlap matrix: must be positive definite
    complex(dp),intent(in)    :: S(:, :)
    !> Refer to [[exp_hermitian_matrix_times_vectors]]
    complex(dp),intent(inout) :: vectors(:, :)
    !> Tolerance for checking if the matrices are hermitian
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

    tolerance = tol_default
    if( present(tol) ) tolerance = tol
    dim = size( H, 1 )
    n_vectors = size( vectors, 2 )
    allocate( eigvecs(dim, n_vectors), aux(dim, n_vectors) )
    allocate( proj(n_vectors, n_vectors) )
    allocate( S_copy, source=S ) 
    allocate( H_copy, source=H )
    allocate( ifail(dim), iwork(5*dim), eigvals(dim), rwork(7*dim) )

    ! Sanity checks
    ! Check if H is hermitian
    call assert( is_hermitian( H, tolerance ), 'H is not hermitian' )
    ! Check if H and vectors have compatible size
    call assert( size( H, 1 ) == size( vectors, 1 ), 'H and vectors have incompatible sizes.' )
    ! Check if S is positive definite
    call assert( is_positive_definite( S, tolerance ), 'S is not positive definite' )
    ! Check if S and vectors have compatible size
    call assert( size( S, 1 ) == size( vectors, 1 ), 'S and vectors have incompatible sizes.' )

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
    call hermitian_matrix_multiply( S, vectors, aux, tol=tolerance )
    call matrix_multiply( eigvecs, aux, proj, 'C')

    ! Now, scale each eigenvector by the exponential of alpha*eigvals
    forall( i = 1:n_vectors ) aux(:, i) = zexp( alpha*eigvals(i) )*eigvecs(:, i)

    call matrix_multiply( aux, proj, vectors )

  end subroutine 


end module matrix_exp