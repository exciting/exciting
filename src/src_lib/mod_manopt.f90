! This file is distributed under the terms of the GNU General Public License.
!
! Considerable parts of this module were strongly adopted from the manopt toolbox
! for matlab developed by Nicolas Boumal and Bamdev Mishra. 
! See http://www.manopt.org
! N. Boumal, B. Mishra, P.-A. Absil, and R. Sepulchre, Journal of Machine Learning Research 15 (2014), 1455-1459
!
! Created October 2019 (Sebastian Tillack)
!
! DESCRIPTION
! This module (together with mod_linalg.f90) provides functionalities for nonlinear optimization
! problems on manifolds. I.e. it provides functions to minimize a real valued cost 
! function f({X}), where {X} is a set of matrices.
! Supported manifolds:
!   real Eulclidean    - X is general real matrix
!   complex Euclidean  - X is general complex matrix
!   real Stiefel       - X is (semi-)orthonormal matrix (X^T*X=1)
!   complex Stiefel    - X is (semi-)unitary matrix (X^H*X=1)
! The individual matrices X do not necessarily have to have the same dimension.
!
! To solve the optimization problem, one of the following routines needs to be invoced:
!
! manopt_MANIFOLD_SOLVER( cost, fun, X, dxo, kx, dx, ...)
!
! where MANIFOLD is
!   euclid  - for Euclidean manifold
!   stiefel - for Stiefel manifold
! and SOLVER is
!   cg      - for a conjugate-gradient solver
!   lbfgs   - for a limited-memory BFGS solver
! The arguments are
!   cost                 - external subroutine to evaluate cost function (matching cost_interface)
!   grad                 - external subroutine to evaluate gradient (matching grad_interface)
!   X(dxo(1),dxo(2),kx)  - on input: initial set of matrices, on exit: minimizer of the cost function
!   dxo(2)               - dimension of the matrices X as allocated in external environment
!   kx                   - number of such X matrices
!   dx(2,kx)             - actual dimension of each X matrix as to use in the optimization
! The following optional arguments (...) are allowed
!   epsgrad     - aimed norm of gradient (default: 1e-6, on exit: aimed norm of gradient)
!   epsfun      - aimed minimal cost function value (default: -infinity, on exit: aimed cost function value)
!   minit       - minimum number of iterations (default: 0, on exit: number of executed iterations)
!   maxit       - maximum number of iterations (default: 1000)
!   minstep     - minimum step length (default: 1e-2)
!   stdout      - output unit for run statistics (default: 0 = no output)
!   eps0        - numerical zero (default: 1e-14)
! for lbfgs solver only
!   memlen      - memory length (default: 12)
!
! Note: The functions manopt_MANIFOLD_SOLVER are overloaded, i.e. their actual implementation
! for real or complex problems is chosen according to the type of X and the interface of cost and grad.
! The user can call it with real or complex quantities and does not need to care.
!
module mod_manopt
  use mod_manopt_matrices
  use mod_manopt_manifolds
  use mod_manopt_solvers

  implicit none

  private

  ! solver interfaces
  interface manopt_euclid_cg
    module procedure :: real_euclid_cg, complex_euclid_cg
  end interface
  interface manopt_euclid_lbfgs
    module procedure :: real_euclid_lbfgs, complex_euclid_lbfgs
  end interface
  interface manopt_stiefel_cg
    module procedure :: real_stiefel_cg, complex_stiefel_cg
  end interface
  interface manopt_stiefel_lbfgs
    module procedure :: real_stiefel_lbfgs, complex_stiefel_lbfgs
  end interface

  ! public solver routines
  public :: manopt_euclid_cg, manopt_euclid_lbfgs, &
            manopt_stiefel_cg, manopt_stiefel_lbfgs

  contains

    !******************************************
    ! USER INTERFACES
    !******************************************
    subroutine real_euclid_cg( X, dxo, kx, dx, &
                               cost, grad, costgrad, update, &
                               epsgrad, epsfun, minit, maxit, stdout, minstep, eps0, noise)
      integer, intent( in)                         :: dxo(2), kx, dx(2,kx)
      real(8), intent( inout)                      :: X(dxo(1),dxo(2),*)
      procedure( real_cost_external), optional     :: cost
      procedure( real_grad_external), optional     :: grad
      procedure( real_costgrad_external), optional :: costgrad
      procedure( real_update_external), optional   :: update
      integer, optional, intent( inout)            :: maxit
      integer, optional, intent( in)               :: minit, stdout
      real(8), optional, intent( inout)            :: epsgrad, epsfun
      real(8), optional, intent( in)               :: minstep, eps0, noise

      type( euclid_manifold) :: M
      type( cg_solver) :: S
      type( real_matrix) :: X0

      ! initialize manifold and X0
      M = euclid_manifold( dxo, kx, dx)
      X0 = real_matrix( X(1:M%DXI(1),1:M%DXI(2),1:M%KX))
      ! initialize solver
      S = cg_solver( M, X0)
      if( present( cost)) S%r_cost_ext => cost
      if( present( grad)) S%r_grad_ext => grad
      if( present( costgrad)) S%r_costgrad_ext => costgrad
      if( present( update)) S%r_update_ext => update
      if( .not. ((associated( S%r_cost_ext) .and. associated( S%r_grad_ext)) .or. associated( S%r_costgrad_ext))) &
        call error( 'real_euclid_cg', 'Either cost and grad or costgrad must be provided.')
      ! set custom parameters
      if( present( epsgrad)) S%epsgrad = epsgrad
      if( present( epsfun))  S%epsfun  = epsfun
      if( present( minit))   S%minit   = minit
      if( present( maxit))   S%maxit   = maxit
      if( present( stdout))  S%stdout  = stdout
      if( present( minstep)) S%minstep = minstep
      if( present( eps0))    S%eps0    = eps0
      if( present( noise))   S%noise   = noise
      ! solve the problem
      call S%solve( M, X0)
      ! assign result
      X(1:M%DXI(1),1:M%DXI(2),1:M%KX) = X0%m
      ! copy further results
      if( present( epsgrad)) epsgrad = S%epsgrad
      if( present( epsfun))  epsfun  = S%cCost
      if( present( maxit))   maxit   = S%maxit
      ! free memory
      call M%clear()
      call S%clear()
      call X0%clear()
    end subroutine real_euclid_cg

    subroutine real_euclid_lbfgs( X, dxo, kx, dx, &
                                  cost, grad, costgrad, update, &
                                  epsgrad, epsfun, minit, maxit, stdout, minstep, eps0, memlen, noise)
      integer, intent( in)                         :: dxo(2), kx, dx(2,kx)
      real(8), intent( inout)                      :: X(dxo(1),dxo(2),*)
      procedure( real_cost_external), optional     :: cost
      procedure( real_grad_external), optional     :: grad
      procedure( real_costgrad_external), optional :: costgrad
      procedure( real_update_external), optional   :: update
      integer, optional, intent( inout)            :: maxit
      integer, optional, intent( in)               :: minit, stdout, memlen
      real(8), optional, intent( inout)            :: epsgrad, epsfun
      real(8), optional, intent( in)               :: minstep, eps0, noise

      type( euclid_manifold) :: M
      type( lbfgs_solver) :: S
      type( real_matrix) :: X0

      ! initialize manifold and X0
      M = euclid_manifold( dxo, kx, dx)
      X0 = real_matrix( X(1:M%DXI(1),1:M%DXI(2),1:M%KX))
      ! initialize solver
      if( present( memlen)) then
        S = lbfgs_solver( M, X0, memlen=memlen)
      else
        S = lbfgs_solver( M, X0)
      end if
      if( present( cost)) S%r_cost_ext => cost
      if( present( grad)) S%r_grad_ext => grad
      if( present( costgrad)) S%r_costgrad_ext => costgrad
      if( present( update)) S%r_update_ext => update
      if( .not. ((associated( S%r_cost_ext) .and. associated( S%r_grad_ext)) .or. associated( S%r_costgrad_ext))) &
        call error( 'real_euclid_lbfgs', 'Either cost and grad or costgrad must be provided.')
      ! set custom parameters
      if( present( epsgrad)) S%epsgrad = epsgrad
      if( present( epsfun))  S%epsfun  = epsfun
      if( present( minit))   S%minit   = minit
      if( present( maxit))   S%maxit   = maxit
      if( present( stdout))  S%stdout  = stdout
      if( present( minstep)) S%minstep = minstep
      if( present( eps0))    S%eps0    = eps0
      if( present( noise))   S%noise   = noise
      ! solve the problem
      call S%solve( M, X0)
      ! assign result
      X(1:M%DXI(1),1:M%DXI(2),1:M%KX) = X0%m
      ! copy further results
      if( present( epsgrad)) epsgrad = S%epsgrad
      if( present( epsfun))  epsfun  = S%cCost
      if( present( maxit))   maxit   = S%maxit
      ! free memory
      call M%clear()
      call S%clear()
      call X0%clear()
    end subroutine real_euclid_lbfgs

    subroutine complex_euclid_cg( X, dxo, kx, dx, &
                                  cost, grad, costgrad, update, &
                                  epsgrad, epsfun, minit, maxit, stdout, minstep, eps0, noise)
      integer, intent( in)                            :: dxo(2), kx, dx(2,kx)
      complex(8), intent( inout)                      :: X(dxo(1),dxo(2),*)
      procedure( complex_cost_external), optional     :: cost
      procedure( complex_grad_external), optional     :: grad
      procedure( complex_costgrad_external), optional :: costgrad
      procedure( complex_update_external), optional   :: update
      integer, optional, intent( inout)               :: maxit
      integer, optional, intent( in)                  :: minit, stdout
      real(8), optional, intent( inout)               :: epsgrad, epsfun
      real(8), optional, intent( in)                  :: minstep, eps0, noise

      type( euclid_manifold) :: M
      type( cg_solver) :: S
      type( complex_matrix) :: X0

      ! initialize manifold and X0
      M = euclid_manifold( dxo, kx, dx)
      X0 = complex_matrix( X(1:M%DXI(1),1:M%DXI(2),1:M%KX))
      ! initialize solver
      S = cg_solver( M, X0)
      if( present( cost)) S%c_cost_ext => cost
      if( present( grad)) S%c_grad_ext => grad
      if( present( costgrad)) S%c_costgrad_ext => costgrad
      if( present( update)) S%c_update_ext => update
      if( .not. ((associated( S%c_cost_ext) .and. associated( S%c_grad_ext)) .or. associated( S%c_costgrad_ext))) &
        call error( 'complex_euclid_cg', 'Either cost and grad or costgrad must be provided.')
      ! set custom parameters
      if( present( epsgrad)) S%epsgrad = epsgrad
      if( present( epsfun))  S%epsfun  = epsfun
      if( present( minit))   S%minit   = minit
      if( present( maxit))   S%maxit   = maxit
      if( present( stdout))  S%stdout  = stdout
      if( present( minstep)) S%minstep = minstep
      if( present( eps0))    S%eps0    = eps0
      if( present( noise))   S%noise   = noise
      ! solve the problem
      call S%solve( M, X0)
      ! assign result
      X(1:M%DXI(1),1:M%DXI(2),1:M%KX) = X0%m
      ! copy further results
      if( present( epsgrad)) epsgrad = S%epsgrad
      if( present( epsfun))  epsfun  = S%cCost
      if( present( maxit))   maxit   = S%maxit
      ! free memory
      call M%clear()
      call S%clear()
      call X0%clear()
    end subroutine complex_euclid_cg

    subroutine complex_euclid_lbfgs( X, dxo, kx, dx, &
                                     cost, grad, costgrad, update, &
                                     epsgrad, epsfun, minit, maxit, stdout, minstep, eps0, memlen, noise)
      integer, intent( in)                           :: dxo(2), kx, dx(2,kx)
      complex(8), intent( inout)                      :: X(dxo(1),dxo(2),*)
      procedure( complex_cost_external), optional     :: cost
      procedure( complex_grad_external), optional     :: grad
      procedure( complex_costgrad_external), optional :: costgrad
      procedure( complex_update_external), optional   :: update
      integer, optional, intent( inout)               :: maxit
      integer, optional, intent( in)                  :: minit, stdout, memlen
      real(8), optional, intent( inout)               :: epsgrad, epsfun
      real(8), optional, intent( in)                  :: minstep, eps0, noise

      type( euclid_manifold) :: M
      type( lbfgs_solver) :: S
      type( complex_matrix) :: X0

      ! initialize manifold and X0
      M = euclid_manifold( dxo, kx, dx)
      X0 = complex_matrix( X(1:M%DXI(1),1:M%DXI(2),1:M%KX))
      ! initialize solver
      if( present( memlen)) then
        S = lbfgs_solver( M, X0, memlen=memlen)
      else
        S = lbfgs_solver( M, X0)
      end if
      if( present( cost)) S%c_cost_ext => cost
      if( present( grad)) S%c_grad_ext => grad
      if( present( costgrad)) S%c_costgrad_ext => costgrad
      if( present( update)) S%c_update_ext => update
      if( .not. ((associated( S%c_cost_ext) .and. associated( S%c_grad_ext)) .or. associated( S%c_costgrad_ext))) &
        call error( 'complex_euclid_lbfgs', 'Either cost and grad or costgrad must be provided.')
      ! set custom parameters
      if( present( epsgrad)) S%epsgrad = epsgrad
      if( present( epsfun))  S%epsfun  = epsfun
      if( present( minit))   S%minit   = minit
      if( present( maxit))   S%maxit   = maxit
      if( present( stdout))  S%stdout  = stdout
      if( present( minstep)) S%minstep = minstep
      if( present( eps0))    S%eps0    = eps0
      if( present( noise))   S%noise   = noise
      ! solve the problem
      call S%solve( M, X0)
      ! assign result
      X(1:M%DXI(1),1:M%DXI(2),1:M%KX) = X0%m
      ! copy further results
      if( present( epsgrad)) epsgrad = S%epsgrad
      if( present( epsfun))  epsfun  = S%cCost
      if( present( maxit))   maxit   = S%maxit
      ! free memory
      call M%clear()
      call S%clear()
      call X0%clear()
    end subroutine complex_euclid_lbfgs

    subroutine real_stiefel_cg( X, dxo, kx, dx, &
                               cost, grad, costgrad, update, &
                               epsgrad, epsfun, minit, maxit, stdout, minstep, eps0, noise)
      integer, intent( in)                         :: dxo(2), kx, dx(2,kx)
      real(8), intent( inout)                      :: X(dxo(1),dxo(2),*)
      procedure( real_cost_external), optional     :: cost
      procedure( real_grad_external), optional     :: grad
      procedure( real_costgrad_external), optional :: costgrad
      procedure( real_update_external), optional   :: update
      integer, optional, intent( inout)            :: maxit
      integer, optional, intent( in)               :: minit, stdout
      real(8), optional, intent( inout)            :: epsgrad, epsfun
      real(8), optional, intent( in)               :: minstep, eps0, noise

      type( stiefel_manifold) :: M
      type( cg_solver) :: S
      type( real_matrix) :: X0

      ! initialize manifold and X0
      M = stiefel_manifold( dxo, kx, dx)
      X0 = real_matrix( X(1:M%DXI(1),1:M%DXI(2),1:M%KX))
      ! initialize solver
      S = cg_solver( M, X0)
      if( present( cost)) S%r_cost_ext => cost
      if( present( grad)) S%r_grad_ext => grad
      if( present( costgrad)) S%r_costgrad_ext => costgrad
      if( present( update)) S%r_update_ext => update
      if( .not. ((associated( S%r_cost_ext) .and. associated( S%r_grad_ext)) .or. associated( S%r_costgrad_ext))) &
        call error( 'real_stiefel_cg', 'Either cost and grad or costgrad must be provided.')
      ! set custom parameters
      if( present( epsgrad)) S%epsgrad = epsgrad
      if( present( epsfun))  S%epsfun  = epsfun
      if( present( minit))   S%minit   = minit
      if( present( maxit))   S%maxit   = maxit
      if( present( stdout))  S%stdout  = stdout
      if( present( minstep)) S%minstep = minstep
      if( present( eps0))    S%eps0    = eps0
      if( present( noise))   S%noise   = noise
      ! solve the problem
      call S%solve( M, X0)
      ! assign result
      X(1:M%DXI(1),1:M%DXI(2),1:M%KX) = X0%m
      ! copy further results
      if( present( epsgrad)) epsgrad = S%epsgrad
      if( present( epsfun))  epsfun  = S%cCost
      if( present( maxit))   maxit   = S%maxit
      ! free memory
      call M%clear()
      call S%clear()
      call X0%clear()
    end subroutine real_stiefel_cg

    subroutine real_stiefel_lbfgs( X, dxo, kx, dx, &
                                  cost, grad, costgrad, update, &
                                  epsgrad, epsfun, minit, maxit, stdout, minstep, eps0, memlen, noise)
      integer, intent( in)                         :: dxo(2), kx, dx(2,kx)
      real(8), intent( inout)                      :: X(dxo(1),dxo(2),*)
      procedure( real_cost_external), optional     :: cost
      procedure( real_grad_external), optional     :: grad
      procedure( real_costgrad_external), optional :: costgrad
      procedure( real_update_external), optional   :: update
      integer, optional, intent( inout)            :: maxit
      integer, optional, intent( in)               :: minit, stdout, memlen
      real(8), optional, intent( inout)            :: epsgrad, epsfun
      real(8), optional, intent( in)               :: minstep, eps0, noise

      type( stiefel_manifold) :: M
      type( lbfgs_solver) :: S
      type( real_matrix) :: X0

      ! initialize manifold and X0
      M = stiefel_manifold( dxo, kx, dx)
      X0 = real_matrix( X(1:M%DXI(1),1:M%DXI(2),1:M%KX))
      ! initialize solver
      if( present( memlen)) then
        S = lbfgs_solver( M, X0, memlen=memlen)
      else
        S = lbfgs_solver( M, X0)
      end if
      if( present( cost)) S%r_cost_ext => cost
      if( present( grad)) S%r_grad_ext => grad
      if( present( costgrad)) S%r_costgrad_ext => costgrad
      if( present( update)) S%r_update_ext => update
      if( .not. ((associated( S%r_cost_ext) .and. associated( S%r_grad_ext)) .or. associated( S%r_costgrad_ext))) &
        call error( 'real_stiefel_lbfgs', 'Either cost and grad or costgrad must be provided.')
      ! set custom parameters
      if( present( epsgrad)) S%epsgrad = epsgrad
      if( present( epsfun))  S%epsfun  = epsfun
      if( present( minit))   S%minit   = minit
      if( present( maxit))   S%maxit   = maxit
      if( present( stdout))  S%stdout  = stdout
      if( present( minstep)) S%minstep = minstep
      if( present( eps0))    S%eps0    = eps0
      if( present( noise))   S%noise   = noise
      ! solve the problem
      call S%solve( M, X0)
      ! assign result
      X(1:M%DXI(1),1:M%DXI(2),1:M%KX) = X0%m
      ! copy further results
      if( present( epsgrad)) epsgrad = S%epsgrad
      if( present( epsfun))  epsfun  = S%cCost
      if( present( maxit))   maxit   = S%maxit
      ! free memory
      call M%clear()
      call S%clear()
      call X0%clear()
    end subroutine real_stiefel_lbfgs

    subroutine complex_stiefel_cg( X, dxo, kx, dx, &
                                  cost, grad, costgrad, update, &
                                  epsgrad, epsfun, minit, maxit, stdout, minstep, eps0, noise)
      integer, intent( in)                            :: dxo(2), kx, dx(2,kx)
      complex(8), intent( inout)                      :: X(dxo(1),dxo(2),*)
      procedure( complex_cost_external), optional     :: cost
      procedure( complex_grad_external), optional     :: grad
      procedure( complex_costgrad_external), optional :: costgrad
      procedure( complex_update_external), optional   :: update
      integer, optional, intent( inout)               :: maxit
      integer, optional, intent( in)                  :: minit, stdout
      real(8), optional, intent( inout)               :: epsgrad, epsfun
      real(8), optional, intent( in)                  :: minstep, eps0, noise

      type( stiefel_manifold) :: M
      type( cg_solver) :: S
      type( complex_matrix) :: X0

      ! initialize manifold and X0
      M = stiefel_manifold( dxo, kx, dx)
      X0 = complex_matrix( X(1:M%DXI(1),1:M%DXI(2),1:M%KX))
      ! initialize solver
      S = cg_solver( M, X0)
      if( present( cost)) S%c_cost_ext => cost
      if( present( grad)) S%c_grad_ext => grad
      if( present( costgrad)) S%c_costgrad_ext => costgrad
      if( present( update)) S%c_update_ext => update
      if( .not. ((associated( S%c_cost_ext) .and. associated( S%c_grad_ext)) .or. associated( S%c_costgrad_ext))) &
        call error( 'complex_stiefel_cg', 'Either cost and grad or costgrad must be provided.')
      ! set custom parameters
      if( present( epsgrad)) S%epsgrad = epsgrad
      if( present( epsfun))  S%epsfun  = epsfun
      if( present( minit))   S%minit   = minit
      if( present( maxit))   S%maxit   = maxit
      if( present( stdout))  S%stdout  = stdout
      if( present( minstep)) S%minstep = minstep
      if( present( eps0))    S%eps0    = eps0
      if( present( noise))   S%noise   = noise
      ! solve the problem
      call S%solve( M, X0)
      ! assign result
      X(1:M%DXI(1),1:M%DXI(2),1:M%KX) = X0%m
      ! copy further results
      if( present( epsgrad)) epsgrad = S%epsgrad
      if( present( epsfun))  epsfun  = S%cCost
      if( present( maxit))   maxit   = S%maxit
      ! free memory
      call M%clear()
      call S%clear()
      call X0%clear()
    end subroutine complex_stiefel_cg

    subroutine complex_stiefel_lbfgs( X, dxo, kx, dx, &
                                     cost, grad, costgrad, update, &
                                     epsgrad, epsfun, minit, maxit, stdout, minstep, eps0, memlen, noise)
      integer, intent( in)                            :: dxo(2), kx, dx(2,kx)
      complex(8), intent( inout)                      :: X(dxo(1),dxo(2),*)
      procedure( complex_cost_external), optional     :: cost
      procedure( complex_grad_external), optional     :: grad
      procedure( complex_costgrad_external), optional :: costgrad
      procedure( complex_update_external), optional   :: update
      integer, optional, intent( inout)               :: maxit
      integer, optional, intent( in)                  :: minit, stdout, memlen
      real(8), optional, intent( inout)               :: epsgrad, epsfun
      real(8), optional, intent( in)                  :: minstep, eps0, noise

      type( stiefel_manifold) :: M
      type( lbfgs_solver) :: S
      type( complex_matrix) :: X0

      ! initialize manifold and X0
      M = stiefel_manifold( dxo, kx, dx)
      X0 = complex_matrix( X(1:M%DXI(1),1:M%DXI(2),1:M%KX))
      ! initialize solver
      if( present( memlen)) then
        S = lbfgs_solver( M, X0, memlen=memlen)
      else
        S = lbfgs_solver( M, X0)
      end if
      if( present( cost)) S%c_cost_ext => cost
      if( present( grad)) S%c_grad_ext => grad
      if( present( costgrad)) S%c_costgrad_ext => costgrad
      if( present( update)) S%c_update_ext => update
      if( .not. ((associated( S%c_cost_ext) .and. associated( S%c_grad_ext)) .or. associated( S%c_costgrad_ext))) &
        call error( 'complex_stiefel_lbfgs', 'Either cost and grad or costgrad must be provided.')
      ! set custom parameters
      if( present( epsgrad)) S%epsgrad = epsgrad
      if( present( epsfun))  S%epsfun  = epsfun
      if( present( minit))   S%minit   = minit
      if( present( maxit))   S%maxit   = maxit
      if( present( stdout))  S%stdout  = stdout
      if( present( minstep)) S%minstep = minstep
      if( present( eps0))    S%eps0    = eps0
      if( present( noise))   S%noise   = noise
      ! solve the problem
      call S%solve( M, X0)
      ! assign result
      X(1:M%DXI(1),1:M%DXI(2),1:M%KX) = X0%m
      ! copy further results
      if( present( epsgrad)) epsgrad = S%epsgrad
      if( present( epsfun))  epsfun  = S%cCost
      if( present( maxit))   maxit   = S%maxit
      ! free memory
      call M%clear()
      call S%clear()
      call X0%clear()
    end subroutine complex_stiefel_lbfgs

    !******************************************
    ! UTILITY FUNCTIONS
    !******************************************
    subroutine error( proc, msg)
      character(*), intent( in) :: proc, msg
      write(*,*)
      write(*,'("Error (mod_manopt/",a,"): ",a)') trim(proc), trim(msg)
      stop
    end subroutine error
end module mod_manopt
