!> Module for solving generalized hermitian eigenvalue problems.
!> The interface combines LAPACK wrappers for
!> DSYGVX and ZHEGVX
module generalized_hermitian_eigenproblem
  use precision, only: dp, i32
  use asserts, only: assert
  use xstring, only: to_char, newline
  use modmpi, only: terminate_if_false
  use math_utils, only: is_hermitian, is_positive_definite
  use lapack_f95_interfaces, only: dsygvx, zhegvx

  implicit none

  private
  public :: solve_generalized_hermitian_eigenproblem

  !> Solve a generalized hermitian-definite eigenproblem
  !> \[
  !>    \mathbf{A} \mathbf{x_i} = \lambda_i \mathbf{B} \mathbf{x_i},
  !> \]
  !> where \( \lambda_i \) is the \(i\)'th eigenvalue and \( \mathbf{x_i} \) the \(i\)'th eigenvector.
  interface solve_generalized_hermitian_eigenproblem
     module procedure :: solve_gen_sym_eigenproblem_real_dp
     module procedure :: solve_gen_hermitian_eigenproblem_complex_dp
  end interface solve_generalized_hermitian_eigenproblem

  !> See [[solve_generalized_hermitian_eigenproblem]].
  !>
  !> This routine acts on the arrays as expected by the LAPACK routines [[dsygvx]] and [[zhegvx]].
  interface xhegvx
     module procedure :: dsygvx_wrapper
     module procedure :: zhegvx_wrapper
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

    num_eigenvalues = size(eigenvalues)

    vl = 0._dp; vu = 0._dp; il = 1

    call xhegvx(1, 'V', 'I', 'U', A, B, vl, vu, il, num_eigenvalues, abstol, eigenvalues_tmp, eigenvectors)
    eigenvalues = eigenvalues_tmp(1:num_eigenvalues)

  end subroutine solve_gen_sym_eigenproblem_real_dp

  !> Same as [[solve_gen_sym_eigenproblem_real_dp]] for complex hermitian matrices
  subroutine solve_gen_hermitian_eigenproblem_complex_dp(A, B, abstol, eigenvalues, eigenvectors)
    !> On entry, the hermitian matrix \( \mathbf{A} \).
    !> On exit, the upper triangle of \( \mathbf{A} \), including the diagonal, is destroyed.
    complex(dp), intent(inout), contiguous :: A(:,:)
    !> On entry, the hermitian, positive definite matrix \( \mathbf{B} \).
    !> On exit, the upper triangle of \( \mathbf{B} \), including the diagonal, is destroyed.
    complex(dp), intent(inout), contiguous :: B(:,:)
    !> See [[solve_gen_sym_eigenproblem_real_dp]]
    real(dp), intent(in) :: abstol
    !> See [[solve_gen_sym_eigenproblem_real_dp]]
    real(dp), intent(out), contiguous:: eigenvalues(:)
    !> Array containing eigenvectors. If not present, only the eigenvalues are calculated
    complex(dp), intent(out), optional, contiguous :: eigenvectors(:,:)

    integer :: il, num_eigenvalues
    real(dp) :: vl, vu
    real(dp), allocatable :: eigenvalues_tmp(:)
    complex(dp) :: eigenvectors_fake(1, 1)

    num_eigenvalues = size(eigenvalues)

    vl = 0._dp; vu = 0._dp; il = 1

    if( present(eigenvectors) ) then
      call xhegvx(1, 'V', 'I', 'U', A, B, vl, vu, il, num_eigenvalues, abstol, eigenvalues_tmp, eigenvectors)
    else
      call xhegvx(1, 'N', 'I', 'U', A, B, vl, vu, il, num_eigenvalues, abstol, eigenvalues_tmp, eigenvectors_fake)
    end if
    eigenvalues = eigenvalues_tmp(1:num_eigenvalues)
  end subroutine

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

    n = size(A, dim=1)
    call assert_char_int_arguments( itype, jobz, range, uplo )
    call assert(is_hermitian(A), 'Matrix A is not hermitian')
    call assert(is_positive_definite(B), 'Matrix B is not positive definite')
    call assert(size(B, dim=1) == n, 'Matrix B dimensions do not match matrix A.')

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

    call terminate_if_diagonalization_failed( info, n, 'dsygvx' )

  end subroutine dsygvx_wrapper

  !> Same as [[dsygvx_wrapper]] for the LAPACK routine [[zhegvx]].
  subroutine zhegvx_wrapper(itype, jobz, range, uplo, A, B, vl, vu, il, iu, abstol, eigenvalues, eigenvectors, lwork_in)
    !> See [[dsygvx_wrapper]]
    integer, intent(in) :: itype
    !> See [[dsygvx_wrapper]]
    character(len=1), intent(in) :: jobz
    !> See [[dsygvx_wrapper]]
    character(len=1), intent(in) :: range
    !> See [[dsygvx_wrapper]]
    character(len=1), intent(in) :: uplo
    !> On entry, the hermitian matrix \( \mathbf{A} \).
    !> On exit, the lower triangle (if `uplo='L'`) or the upper triangle (if `uplo='U'`) of \( \mathbf{A} \),
    !> including the diagonal, is destroyed.
    complex(dp), intent(inout), contiguous :: A(:,:)
    !> On entry, the positive definite matrix \( \mathbf{B} \).
    !> On exit, the lower triangle (if `uplo='L'`) or the upper triangle (if `uplo='U'`) of \( \mathbf{B} \),
    !> including the diagonal, is destroyed.
    complex(dp), intent(inout), contiguous :: B(:,:)
    !> See [[dsygvx_wrapper]]
    real(dp), intent(in):: vl
    !> See [[dsygvx_wrapper]]
    real(dp), intent(in):: vu
    !> See [[dsygvx_wrapper]]
    integer(i32), intent(in) :: il
    !> See [[dsygvx_wrapper]]
    integer(i32), intent(in) :: iu
    !> See [[dsygvx_wrapper]]
    real(dp), intent(in) :: abstol
    !> See [[dsygvx_wrapper]]
    real(dp), intent(out), allocatable :: eigenvalues(:)
    !> Array containing eigenvectors. 
    complex(dp), intent(out), contiguous :: eigenvectors(:,:)
    !> See [[dsygvx_wrapper]]
    integer(i32), intent(in), optional :: lwork_in

    integer(i32) :: m, lwork, info, n, ldz, i
    real(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: work(:)
    integer(i32), allocatable :: ifail(:), iwork(:)

    n = size(A, dim=1)
    call assert_char_int_arguments( itype, jobz, range, uplo )
    call assert(is_hermitian(A), 'Matrix A is not hermitian')
    call assert(is_positive_definite(B), 'Matrix B is not positive definite')
    call assert(size(B, dim=1) == n, 'Matrix B dimensions do not match matrix A.')

    ldz = size(eigenvectors, dim=1)

    allocate(iwork(5*n))
    allocate(ifail(n))
    allocate(eigenvalues(n))
    allocate(rwork(7*n))

    if (present(lwork_in)) then
       lwork = lwork_in
    else
       lwork = -1; allocate(work(1))
       call zhegvx(itype, jobz, range, uplo, n, A, n, B, n, vl, vu, il, iu, abstol, m, eigenvalues, eigenvectors, ldz, &
            work, lwork, rwork, iwork, ifail, info)
       call terminate_if_false(info == 0, 'dsygvx work space query failed.')
       lwork = work(1); deallocate(work)
    end if

    allocate(work(lwork))

    call zhegvx(itype, jobz, range, uplo, n, A, n, B, n, vl, vu, il, iu, abstol, m, eigenvalues, eigenvectors, ldz, &
         work, lwork, rwork, iwork, ifail, info)

    call terminate_if_diagonalization_failed( info, n, 'zhegvx' )

  end subroutine

  !> (private) Assertions required by `xhegvx` on the character and integers arguments
  subroutine assert_char_int_arguments( itype, jobz, range, uplo )
    !> See description in [[dsygvx_wrapper]]
    integer, intent(in) :: itype
    !> See description in [[dsygvx_wrapper]]
    character(len=1), intent(in) :: jobz
    !> See description in [[dsygvx_wrapper]]
    character(len=1), intent(in) :: range
    !> See description in [[dsygvx_wrapper]]
    character(len=1), intent(in) :: uplo

    call assert(any(itype == allowed_itype), 'itype is not one of the allowed characters (1, 2, 3).')
    call assert(any(jobz == allowed_jobz), 'jobz is not one of the allowed characters ("N", "V").')
    call assert(any(range == allowed_range), 'range is not one of the allowed characters ("A", "V", "I").')
    call assert(any(uplo == allowed_uplo), 'uplo is not one of the allowed characters ("U", "L").')
  end subroutine

  !> (private) Terminate if `info` is not equal to zero (which means that the diagonalization failed)
  subroutine terminate_if_diagonalization_failed( info, n, xhegvx_name )
    !> Variable that decides if the code must be terminated
    integer(i32), intent(in) :: info
    !> Matrix dimension
    integer(i32), intent(in) :: n
    !> Diagonalization subroutine that failed
    character(len=*), intent(in) :: xhegvx_name

    character(:), allocatable :: error_message

    if( info /= 0 ) then
      error_message = xhegvx_name // ' routine failed with INFO = ' // to_char(info)
      if (info > n) then
        error_message = error_message // newline // 'The leading minor of the overlap matrix of order ' // &
              to_char( info - n ) // ' is not positive definite.' // newline // &
              'Order of overlap matrix: ' // to_char(n)
      end if
      call terminate_if_false( .false., error_message )
    end if
  end subroutine

end module generalized_hermitian_eigenproblem
