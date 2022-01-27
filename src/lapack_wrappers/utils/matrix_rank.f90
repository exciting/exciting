!> Module for calculating the matrix rank with
module matrix_rank
  use precision, only: dp
  use math_utils, only: all_zero
  use singular_value_decomposition, only: svd_divide_conquer

  implicit none

  private
  public :: matrix_rank_SVD


  !> Tolerance for defining \(\sigma_i = 0 \)
  real(dp), parameter :: default_tol_sigma = 1e-10_dp


  !> Determine the rank of a matrix \( \mathbf{A} \) by calculating the SVD and
  !> counting the non-zero singular values.
  !> A singular value is understood as non-zero if it is
  !> greater than a tolerance **tol**.
  interface matrix_rank_SVD
    module procedure matrix_rank_SVD_real_dp, matrix_rank_SVD_complex_dp
  end interface matrix_rank_SVD

  contains
  !---------------------------------------------------------------------------------------------------------------------
  ! matrix_rank_SVD

  !> Calculate the rank of a matrix \( \mathbf{A} \) by determining
  !> the singular value decomposition. The rank is
  !> given by the number of non-zero singular values.
  !> A singular value is understood as non-zero if it is
  !> greater than a tolerance **tol**.
  function matrix_rank_SVD_real_dp(A, tol) result(rank)
    !> Input matrix \( \mathbf{A} \)
    real(dp), intent(in), contiguous :: A(:, :)
    !> Tolerance for defining
    !> \(\sigma_i \equiv 0 \Leftrightarrow \sigma_i \le \text{tol} \)
    real(dp), intent(in), optional :: tol

    integer :: rank

    real(dp) :: tol_
    real(dp), allocatable :: sigma(:)
    integer :: i

    tol_ = default_tol_sigma
    if(present(tol)) tol_ = tol

    allocate(sigma(min(size(A, dim=1), size(A, dim=2))))
    call svd_divide_conquer(A, sigma)

    rank = size(sigma)
    do i=size(sigma), 1, -1
      if (all_zero(sigma(i), tol=tol_)) then
        rank = rank - 1
      else
        exit
      end if
    end do
  end function matrix_rank_SVD_real_dp


  !> Calculate the rank of a matrix \( \mathbf{A} \) by determining
  !> the singular value decomposition. The rank is
  !> given by the number of non-zero singular values.
  !> A singular value is understood as non-zero if it is
  !> greater than a tolerance **tol**.
  function matrix_rank_SVD_complex_dp(A, tol) result(rank)
    !> Input matrix \( \mathbf{A} \)
    complex(dp), intent(in), contiguous :: A(:, :)
    !> Tolerance for defining
    !> \(\sigma_i \equiv 0 \Leftrightarrow \sigma_i \le \text{tol} \)
    real(dp), intent(in), optional :: tol

    integer :: rank

    real(dp) :: tol_
    real(dp), allocatable :: sigma(:)
    integer :: i

    tol_ = default_tol_sigma
    if(present(tol)) tol_ = tol

    allocate(sigma(min(size(A, dim=1), size(A, dim=2))))

    call svd_divide_conquer(A, sigma)

    rank = size(sigma)
    do i=size(sigma), 1, -1
      if (all_zero(sigma(i), tol=tol_)) then
        rank = rank - 1
      else
        exit
      end if
    end do
  end function matrix_rank_SVD_complex_dp

end module matrix_rank