!> Lanczos algorithm for the ISDF BSH as implemented for fast BSE.
module iterative_solver
  use precision, only: dp, i32
  use asserts, only: assert
  use modmpi, only: terminate_if_false
  use math_utils, only: all_zero
  use xlapack, only: norm

  implicit none

  private
  public :: lanczos

  !> Zero tolerance
  real(dp), parameter :: zero_tol = 1e-10

  contains

  !> Calculate \( k \) lanczos steps for the Bethe-Salpeter Hamiltonian \( \mathbf{H} \) such that
  !> \[
  !>   \mathbf{H} \cdot \mathbf{Q}_k = \mathbf{Q}_k \cdot \mathbf{T}_k
  !> \]
  !> where \( \mathbf{Q}_k^H \mathbf{Q}_k = \mathbf I_{k \times k}, \: \mathbf{Q} \in \mathcal{C}^{N \times k} \) is,
  !> a transformation matrix, \( \mathbf{T}_k \in \mathcal{R}^{k \times \k} \) is a symmetric tridiagonal matrix,
  !> and \( N \) is the size of the Hamiltonian. \( k \) is chosen such that \( k \ll N \).
  !> The diagonal parts \( t_{ii} \) and sub diagonal parts \( t_{ii+1} \) are calculated iteratively as
  !> \[
  !>   t_{ii} = \mathbf{q_i}^H \cdot \mathbf{H} \mathbf{q_i},
  !> \]
  !> \[
  !>   t_{i+1i} = || (\mathbf{H} - t_{ii}) \cdot \mathbf{q_i} ||,
  !> \]
  !> where \( \mathbf{q_i} \) is the \( i \)'th column of \( \mathbf{Q} \). \( \mathbf{q}_{i+1} \) is calculated 
  !> as follows
  !> \[
  !>   \mathbf{q}_{i+1} = \frac{ \mathbf{H} - t_{ii}}{t_{i+1}} \cdot \mathbf{q_i}.
  !> \]
  !> The first column of \( \mathbf{Q} \), \( \mathbf{q}_1 \) must be given.
  !> If the algorithm breaks down before the \( k \)'th itereration it returns the results so far.
  !> If the algorithm breaks down in the first iteration, the output arrays stay deallocated.
  !> It is the responsibility of the user to verify that [[lanczos]] did not fail.
  subroutine lanczos(k, matrix_vector_product, q_1, alpha, beta, Q_k)
    !> Maximum number of lanczos iterations
    integer, intent(in) :: k
    !> Matrix vector product to be used
    interface
      subroutine matrix_vector_product(vector_in, vector_out)
        use precision, only: dp
        complex(dp), intent(in) :: vector_in(:)
        complex(dp), intent(out) :: vector_out(:)
      end subroutine 
    end interface 
    !> First column of \( \mathbf{Q}_k \), needs to be normalized.
    complex(dp), intent(in) :: q_1(:)
    !> Diagonal of \( \mathbf{T}_k \)
    real(dp), intent(out), allocatable :: alpha(:)
    !> Sub diagonal of \( \mathbf{T}_k \)
    real(dp), intent(out), allocatable :: beta(:)
    !> Transformation matrix. Will be only saved if present.
    complex(dp), intent(out), allocatable, optional :: Q_k(:, :)
    
    logical :: save_Q
    integer(i32) :: iter, k_run, n_matrix
    real(dp), allocatable :: alpha_(:), beta_(:)
    complex(dp), allocatable :: x(:), q_vec(:), q_vec_old(:)

    n_matrix = size(q_1)
    save_Q = present(Q_k)

    call assert(k >= 1, 'k is smaller then 1..')
    call assert(k <= n_matrix, 'k is larger than the size of the matrix.')
    call assert(norm(q_1) >= zero_tol, 'norm(q_1) is zero.')

    allocate(alpha(k))
    allocate(beta(0 : k))
    allocate(x(n_matrix))
    if(save_Q) allocate(Q_k(n_matrix, k+1))

    q_vec = q_1 / norm(q_1)
    q_vec_old = q_vec
    if (save_Q) Q_k(:, 1) = q_vec

    beta(0) = 0.0_dp
    k_run = 0
    do iter=1, k
      call matrix_vector_product(q_vec, x)
      x = x - beta(iter-1) * q_vec_old
      alpha(iter) = real(dot_product(q_vec, x), kind=dp)
      x = x - alpha(iter) * q_vec
      beta(iter) = norm(x)

      ! Break loop if linear independence is reached
      if (beta(iter) <= zero_tol) exit

      q_vec_old = q_vec
      q_vec = x / beta(iter)
      if (save_Q) Q_k(:, iter+1) = q_vec
      k_run = iter
    end do

    if(k_run > 0) then
      alpha = alpha(:k_run)
      beta = beta(1 : k_run)
      if(save_Q) Q_k = Q_k(:, :k_run)
    else
      ! Break down in the first iteration cannot yield a result
      deallocate(alpha, beta)
      if(save_Q) deallocate(Q_k)
    end if
  end subroutine lanczos

end module iterative_solver