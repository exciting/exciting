module integration
  use asserts, only: assert
  use constants, only: zone, zzero
  use hermitian_matrix_multiplication, only: hermitian_matrix_multiply
  use math_utils, only: is_hermitian, is_positive_definite
  use precision, only: dp

  implicit none
  private
  public :: ODESolver_RungeKutta4thOrder

  contains 

  !> Runge-Kutta of 4th order to calculate \( x(t+\Delta t) \)
  !> following
  !>\[
  !>	   \alpha S\frac{d}{dt}| x(t)\rangle = H(t)| x(t) \rangle,
  !>  \]
  !> provided \( x(t) \).
  !> The basic formulas employed here are
  !> \[
  !>	   | x(t+\Delta t)\rangle = | x(t) \rangle
  !>    	+ \frac{\Delta t}{6 \alpha}(k_1+2k_2+2k_3+k_4)
  !>  \]
  !>	where
  !> \[
  !>	   k_1 = S^{-1}H(t)| x(t) \rangle,
  !>  \]
  !> \[
  !>     k_2 = S^{-1}H\left(t + \frac{\Delta t}{2}\right)
  !>    \left[| x(t) \rangle + k_1 \frac{\Delta t}{2}\right],
  !>  \]
  !> \[
  !>	   k_3 = S^{-1}H\left(t + \frac{\Delta t}{2}\right)
  !>    \left[| x(t) \rangle + k_2 \frac{\Delta t}{2}\right],
  !>  \]
  !> \[
  !>     k_4 = S^{-1}H(t+\Delta t)\left[| x(t) \rangle + k_3\Delta t \right].
  !>  \]
  subroutine ODESolver_RungeKutta4thOrder( time_step, alpha, &
    & H, H_past, S, x )
    !> Time step
    real(dp), intent(in)          :: time_step
    !> Complex prefactor
    complex(dp), intent(in)       :: alpha
    !> Hamiltonian matrix at time \( t \) 
    complex(dp), intent(in)       :: H(:, :)
    !> Hamiltonian matrix at time \( t - \Delta t \)  
    complex(dp), intent(in)       :: H_past(:, :)
    !> Overlap matrix
    complex(dp), intent(in)       :: S(:, :)
    !> Vectors to evolve
    complex(dp), intent(inout)    :: x(:, :)
  
    integer                       :: i, info
    integer                       :: dim, n_vectors
    complex(dp)                   :: prefactor
    complex(dp), allocatable      :: k(:, :, :), y(:, :), H_aux(:, :), S_aux(:, :)
  
  
    dim = size( H, 1 )
    n_vectors = size( x, 2 )
    ! Sanity checks
    ! Check if H and H_past are hermitian
    call assert( is_hermitian( H ), 'H is not hermitian' )
    call assert( is_hermitian( H_past ), 'H_past is not hermitian' )
    ! Check if H, H_past and x have compatible size
    call assert( size( H, 1 ) == size( x, 1 ), 'H and x have incompatible sizes.' )
    call assert( size( H_past, 1 ) == size( x, 1 ), 'H_past and x have incompatible sizes.' )
    ! Check if S is positive definite
    call assert( is_positive_definite( S ), 'S is not positive definite' )
    ! Check if S and x have compatible size
    call assert( size( S, 1 ) == size( x, 1 ), 'S and x have incompatible sizes.' )
  
    ! Allocate  
    allocate( k(dim, n_vectors, 4), y(dim, n_vectors) )
    allocate( H_aux(dim, dim), S_aux(dim, dim) )
  
    ! Initiliaze
    prefactor = time_step/alpha
    y = x
    S_aux = S
    H_aux = H
  
    ! Obtain k_1, ..., k_4
    do i = 1, 4
      ! k(:, :, i) = prefactor*H_aux*y, H_aux must be hermitian
      call hermitian_matrix_multiply( H_aux, prefactor*y, k(:, :, i) )
      ! Obtain (S^(-1))*k(:, :, i) for positive definite S (k will store the solution)
      call ZPOSV( 'U', dim, n_vectors, S_aux, dim, k(:, :, i), dim, info )
      ! Restores S_aux to its original value, after being modified by ZPOSV
      S_aux = S
  
      select case( i )
        case( 1 )
          y = x + 0.5_dp*k(:, :, 1)
          ! H_aux = an estimate of H(t+0.5*\Delta t)
          H_aux = 1.5_dp*H - 0.5_dp*H_past
        case( 2 )
          y = x + 0.5_dp*k(:, :, 2)
        case( 3 )
          y = x + k(:, :, 3)
          ! H_aux = an estimate of H(t+\Delta t)
          H_aux = 2._dp*H -  H_past
        case( 4 )
          exit
      end select
    end do
  
    x = x + (1._dp/6)*( k(:, :, 1) + 2._dp*k(:, :, 2) + 2._dp*k(:, :, 3) + k(:, :, 4) )
  
  end subroutine ODESolver_RungeKutta4thOrder

end module