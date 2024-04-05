!> Lanczos algorithm for the ISDF BSH as implemented for fast BSE.
module iterative_solver
  use precision, only: dp
  use asserts, only: assert
  use math_utils, only: all_zero
  use xlapack, only: norm

  implicit none

  private
  public :: lanczos

  contains

  !> Calculate \( k \) lanczos steps for the Bethe-Salpeter Hamiltonian \( \mathbf{H} \) such that
  !> \[
  !>   \mathbf{H} \cdot \mathbf{Q}_k = \mathbf{Q}_k \cdot \mathbf{T}_k
  !> \]
  !> where \( \mathbf{Q}_k^H \mathbf{Q}_k = \mathbf I_{k \times k}, \mathbf{Q} \in \mathcal{C}^{N \times k \) is,
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
  !> [
  !>   \mathbf{q}_{i+1} = \frac{ \mathbf{H} - t_{ii}}{t_{i+1}} \cdot \mathbf{q_i}.
  !> ]
  !> The first column of \( \mathbf{Q} \), \( \mathbf{q}_1 \) must be given.
  !> If the algorithm breaks down before the \( k \)'th itereration it returns the results so far.
  subroutine lanczos(k, matrix_vector_product, q_1, alpha, beta, save_Q, Q_k)
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
    !> Save the transformation matrix. If set to false, `[[Q_k]]` is not
    !> touched.
    logical, intent(in) :: save_Q
    !> Transformation matrix. If `[[save_Q]]` is `.false.`, this matrix will not be allocated.
    complex(dp), intent(out), allocatable :: Q_k(:, :)
    
    ! local variables
    integer :: iter, k_, n_matrix 
    real(dp), allocatable :: alpha_(:), beta_(:)
    complex(dp), allocatable :: Q_k_(:, :), q_vec(:), q_vec_old(:), x(:)

    n_matrix = size(q_1)
    call assert(k <= n_matrix, 'k is larger than the size of the matrix.')

    allocate(alpha_(k))
    allocate(beta_(0 : k))
    allocate(x(n_matrix))

    q_vec = q_1 / norm(q_1)
    q_vec_old = q_vec
    beta_(0) = 0.0_dp

    if (save_Q) then
      allocate(Q_k_(n_matrix, k+1))
      Q_k_(:, 1) = q_vec
    end if

    
    do iter=1, k
      call matrix_vector_product(q_vec, x)
      x = x - beta_(iter-1) * q_vec_old
      alpha_(iter) = real(dot_product(q_vec, x), kind=dp)
      x = x - alpha_(iter) * q_vec
      beta_(iter) = norm(x)

      ! Break loop if linear independence is reached
      if (all_zero(beta_(iter))) exit

      q_vec_old = q_vec
      q_vec = x / beta_(iter)
      if (save_Q) Q_k_(:, iter+1) = q_vec
      k_ = iter
    end do

    alpha = alpha_(:k_)
    beta = beta_(1 : k_)
    if (save_Q) Q_k = Q_k_(:, :k_)
  end subroutine lanczos

end module iterative_solver