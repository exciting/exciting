!> Module for solving generalized symmetric definite eigenvalue problems.
!> The interface combines LAPACK wrappers for
!> DSYGVX
module generalized_symmetric_eigenproblem
  use precision, only: dp
  use asserts, only: assert
  use xstring, only: to_char, newline
  use modmpi, only: terminate_if_false
  use math_utils, only: is_hermitian
  use lapack_f95_interfaces, only: dsygvx

  implicit none

  private
  public :: solve_generalized_symmetric_eigenproblem, xhegvx

  !> Solve a generalized symmetric-definite eigenproblem
  !> \[
  !>    \mathbf{A} \mathbf{x_i} = \lambda_i \mathbf{B} \mathbf{x_i},
  !> \]
  !> where \( \lambda_i \) is the \(i\)'th eigenvalue and \( \mathbf{x_i} \) the \(i\)'th eigenvector.
  interface solve_generalized_symmetric_eigenproblem
     module procedure :: solve_gen_sym_eigenproblem_real_dp
  end interface solve_generalized_symmetric_eigenproblem

  !> See [[solve_generalized_symmetric_eigenproblem]].
  !>
  !> This routine acts on the arrays as expected by the LAPACK routine [[dsygvx]].
  interface xhegvx
     module procedure dsygvx_wrapper
  end interface xhegvx

  !> Allowed input for `itype`
  integer, parameter :: allowed_itype(3) = [1, 2, 3]
  !> Allowed input for `jobz`
  character(1), parameter :: allowed_jobz(2) = ['N', 'V']
  !> Allowed input for `range`
  character(1), parameter :: allowed_range(3) = ['A', 'V', 'I']
  !> Allowed input for `uplo`
  character(1), parameter :: allowed_uplo(2) = ['U', 'L']

contains

  !> Solve a real symmetric-definite generalized eigenvalue problem of the form
  !> \[
  !>    \mathbf{A} \mathbf{x_i} = \lambda_i \mathbf{B} \mathbf{x_i},
  !> \]
  !> where \( \mathbf{A} \) and \( \mathbf{B} \) are real-symmetric matrices, and \( \mathbf{B} \)
  !> is positive definite.
  !> This routine computes eigenvalues as well as eigenvectors.
  !> The number of computed eigenvalues corresponds to the size of the array provided in the `eigenvalues` argument.
  subroutine solve_gen_sym_eigenproblem_real_dp(A, B, abstol, eigenvalues, eigenvectors)
    !> On entry, the real-symmetric matrix \( \mathbf{A} \).
    !> On exit, the upper triangle of \( \mathbf{A} \), including the diagonal, is destroyed.
    real(dp), intent(inout), contiguous :: A(:,:)
    !> On entry, the real-symmetric, positive definite matrix \( \mathbf{B} \).
    !> On exit, the upper triangle of \( \mathbf{B} \), including the diagonal, is destroyed.
    real(dp), intent(inout), contiguous :: B(:,:)
    !> Absolute error tolerance for eigenvalues.
    real(dp), intent(in) :: abstol
    !> Array containing eigenvalues in ascending order.
    real(dp), intent(out), contiguous:: eigenvalues(:)
    !> Array containing eigenvectors.
    real(dp), intent(out), contiguous :: eigenvectors(:,:)

    integer :: il, num_eigenvalues
    real(dp) :: vl, vu
    real(dp), allocatable :: eigenvalues_tmp(:)

    call assert(is_hermitian(A), 'Matrix A is not symmetric.')
    call assert(is_hermitian(B), 'Matrix B is not symmetric.')

    num_eigenvalues = size(eigenvalues)

    vl = 0._dp; vu = 0._dp; il = 1

    call xhegvx(1, 'V', 'I', 'U', A, B, vl, vu, il, num_eigenvalues, abstol, eigenvalues_tmp, eigenvectors)
    eigenvalues = eigenvalues_tmp(1:num_eigenvalues)


  end subroutine solve_gen_sym_eigenproblem_real_dp

  !> Solves a real generalized symmetric-definite eigenvalue problem
  !> where \( \mathbf{A} \) and \( \mathbf{B} \) are \( N \times N \) symmetric matrices, and \( \mathbf{B} \)
  !> is positive definite.
  !> This routine acts on the arrays as expected by the LAPACK routine [[dsygvx]].
  subroutine dsygvx_wrapper(itype, jobz, range, uplo, A, B, vl, vu, il, iu, abstol, eigenvalues, eigenvectors, lwork_in)
    !> Specifies the problem type to be solved by DSYGVX:
    !> `1`: \( \mathbf{A} \mathbf{x} = \lambda \mathbf{B} \mathbf{x} \),
    !> `2`: \( \mathbf{A} \mathbf{B} \mathbf{x} = \lambda \mathbf{x} \),
    !> `3`: \( \mathbf{B} \mathbf{A} \mathbf{x} = \lambda \mathbf{x} \).
    integer, intent(in) :: itype
    !> Specifies whether to compute eigenvectors:
    !> `'N'`: Compute eigenvalues only.
    !> `'V'`: Compute eigenvalues and eigenvectors.
    character(len=1), intent(in) :: jobz
    !> Specifies the range of eigenvalues to be found:
    !> `'A'`: All eigenvalues will be found.
    !> `'V'`: Eigenvalues in the interval \( (vl, vu] \) will be found.
    !> `'I'`: Eigenvalues indexed from `il` to `iu` will be found.
    character(len=1), intent(in) :: range
    !> Specifies whether the upper or lower triangular part of matrices is stored:
    !> `'U'`: Upper triangular part is stored.
    !> `'L'`: Lower triangular part is stored.
    character(len=1), intent(in) :: uplo
    !> On entry, the symmetric matrix \( \mathbf{A} \).
    !> On exit, the lower triangle (if `uplo='L'`) or the upper triangle (if `uplo='U'`) of \( \mathbf{A} \),
    !> including the diagonal, is destroyed.
    real(dp), intent(inout), contiguous :: A(:,:)
    !> On entry, the symmetric matrix \( \mathbf{B} \).
    !> On exit, the lower triangle (if `uplo='L'`) or the upper triangle (if `uplo='U'`) of \( \mathbf{B} \),
    !> including the diagonal, is destroyed.
    real(dp), intent(inout), contiguous :: B(:,:)
    !> If `range = 'V'`, specifies the lower bound of the interval to search for eigenvalues.
    real(dp), intent(in):: vl
    !> If `range = 'V'`, specifies the upper bound of the interval to search for eigenvalues.
    real(dp), intent(in):: vu
    !> If `range = 'I'`, specifies the index of the smallest eigenvalue to be returned.
    integer, intent(in) :: il
    !> If `range = 'I'`, specifies the index of the smallest eigenvalue to be returned.
    integer, intent(in) :: iu
    !> Absolute error tolerance for eigenvalues.
    real(dp), intent(in) :: abstol
    !> Array containing eigenvalues in ascending order.
    real(dp), intent(out), allocatable :: eigenvalues(:)
    !> Array containing eigenvectors.
    real(dp), intent(out), contiguous :: eigenvectors(:,:)
    !> Dimension of the work array. If not present, a work space query is done.
    integer, intent(in), optional :: lwork_in

    integer :: m, lwork, info, n, ldz, i
    real(dp), allocatable :: work(:)
    integer, allocatable :: ifail(:), iwork(:)
    character(:), allocatable :: error_message

    call assert(any(itype == allowed_itype), 'itype is not one of the allowed characters (1, 2, 3).')
    call assert(any(jobz == allowed_jobz), 'jobz is not one of the allowed characters ("N", "V").')
    call assert(any(range == allowed_range), 'range is not one of the allowed characters ("A", "V", "I").')
    call assert(any(uplo == allowed_uplo), 'uplo is not one of the allowed characters ("U", "L").')

    n = size(A, dim=1)
    call assert(size(A, dim=2) == n, 'Matrix A is not square.')
    call assert(size(B, dim=1) == n .and. size(B, dim=2) == n, 'Matrix B dimensions do not match matrix A.')

    ldz = size(eigenvectors, dim=1)

    allocate(iwork(5*n))
    allocate(ifail(n))
    allocate(eigenvalues(n))

    if (present(lwork_in)) then
       lwork = lwork_in
    else
       lwork = -1; allocate(work(1))
       call dsygvx(itype, jobz, range, uplo, n, A, n, B, n, vl, vu, il, iu, abstol, m, eigenvalues, eigenvectors, ldz, &
            work, lwork, iwork, ifail, info)
       call terminate_if_false(info == 0, 'dsygvx work space query failed.')
       lwork = work(1); deallocate(work)
    end if

    allocate(work(lwork))

    call dsygvx(itype, jobz, range, uplo, n, A, n, B, n, vl, vu, il, iu, abstol, m, eigenvalues, eigenvectors, ldz, &
         work, lwork, iwork, ifail, info)

    error_message = 'dsygvx routine failed with INFO = ' // to_char(info)

    if (info > n) then
       i = info - n
       error_message = error_message // newline // 'The leading minor of the overlap matrix of order ' // &
            to_char(i) // ' is not positive definite.' // newline // &
            'Order of overlap matrix: ' // to_char(n)
    end if

    call terminate_if_false(info == 0, error_message)

  end subroutine dsygvx_wrapper

end module generalized_symmetric_eigenproblem
