!> Module for diagonalizing tridiagonal matrices.
!> The interface combines LAPACK wrappers for 
!> DSTEDC
module diagonalize_tridiagonal
  use precision, only: dp
  use asserts, only: assert
  use modmpi, only: terminate_if_false
  use math_utils, only: is_square, is_unitary
  use lapack_f95_interfaces, only: dstedc

  implicit none

  private
  public :: diagonalize_symtridiag, xstedc

  !> Diagonalize a tridiagonal symmetric matrix \( \mathbf{T} \) such that
  !> \[
  !>    \mathbf{T} \mathbf{x_i} = \lambda_i \mathbf{x_i},
  !> \]
  !> where \( \lambda_i \) is the \(i\)'th eigenvalue and \( \mathbf{x_i} \) \(i\)'th eigenvector.
  interface diagonalize_symtridiag
    module procedure :: diagonalize_symtridiag_real_dp
  end interface diagonalize_symtridiag

  !> See [[diagonalize_symtridiag]].
  !>
  !> This routine acts on the arrays as expected by the LAPACK routine [[dstedc]].
  interface xstedc
    module procedure dstedc_wrapper
  end interface xstedc

  !> See [[diagonalize_symtridiag]].
  !>
  !> This routine acts on the arrays as expected by the LAPACK routine [[dstedc]].
  interface xstedc_wsq
    module procedure dstedc_wrapper_workspace_query
  end interface xstedc_wsq

  !> Allowed input for `compz`
  character(1), parameter :: allowed_compz(3) = ['N', 'I', 'V']


  contains
  
  !> Diagonalize a real tridiagonal symmetric matrix \( \mathbf{T} \) such that
  !> \[
  !>    \mathbf{T} \mathbf{x_i} = \lambda_i \mathbf{x_i},
  !> \]
  !> where \( \lambda_i \) is the \(i\)'th eigenvalue and \( \mathbf{x_i} \) \(i\)'th eigenvector.
  subroutine diagonalize_symtridiag_real_dp( &
    diagonal_in, &
    subdiagonal_in, &
    eigenvalues_out, &
    eigenvectors_out, &
    transformation_matrix_in &
  )
    !> Diagonal of \( \mathbf{T} \)
    real(dp), intent(in) :: diagonal_in(:)
    !> Subdiagonal of \( \mathbf{T} \)
    real(dp), intent(in) :: subdiagonal_in(:)
    !> On exit: Eigenvalues of \( \mathbf{T} \)
    real(dp), intent(out), allocatable :: eigenvalues_out(:)
    !> If present on exit: Eigenvectors of \( \mathbf{T} \), columb-wise stored.
    !> If [[transformation_matrix]] is given, contains the eigenvectors of the symmetric matrix  \(\mathbf{M}\),
    !> that was transformed to \(\mathbf{T}\) (see [[tranformation_matrix]] for detailed information on this option).
    real(dp), intent(out), allocatable, optional :: eigenvectors_out(:, :)
    !> Orthogonal matrix \(\mathbf{Q}\) used to transform the symmetric matrix \(\mathbf{M}\) to the tridiagonal matrix
    !> \( \mathbf{T} \) such that:
    !> \[
    !>   \mathbf{T} = \mathbf{Q^\top} \mathbf{M} \mathbf{Q}.
    !> \]
    !> If present and `[[eigenvectors]]` is present, `eigenvectors` contains on output the eigenvectors of the
    !> symmetric matrix. If `eigenvectors` is not present, this matrix will be ignored if given.
    real(dp), intent(in), optional :: transformation_matrix_in(:, :)

    integer :: N
    real(dp), allocatable :: subdiagonal(:), eigenvectors_dummy(:, :)

    N = size(diagonal_in)
    call assert(size(subdiagonal_in) == N-1, 'Size of subdiagonal is not smaller by 1 then size of diagonal.')
    eigenvalues_out = diagonal_in
    subdiagonal = subdiagonal_in

    ! Calculate only eigenvalues
    if (.not. present(eigenvectors_out)) then
      allocate(eigenvectors_dummy(N, N))
      call xstedc('N', N, eigenvalues_out, subdiagonal, eigenvectors_dummy)

    ! Calculate eigenvalues and eigenvalues of the tridiagonal matrix
    elseif(present(eigenvectors_out) .and. .not. present(transformation_matrix_in)) then
      if (allocated(eigenvectors_out)) deallocate(eigenvectors_out); allocate(eigenvectors_out(N, N))
      call xstedc('I', N, eigenvalues_out, subdiagonal, eigenvectors_out)

    ! Calculate eigenvalues and eigenvectors of the symmetric matrix that was tridiagonalized
    elseif(present(eigenvectors_out) .and. present(transformation_matrix_in)) then
      call assert(all(shape(transformation_matrix_in) == [N, N]), 'transformation_matrix has not shape [N, N].')
      call assert(is_unitary(transformation_matrix_in), 'tranformation_matrix is not orthogonal.')
      eigenvectors_out = transformation_matrix_in
      call xstedc('V', N, eigenvalues_out, subdiagonal, eigenvectors_out)
    end if
  end subroutine diagonalize_symtridiag_real_dp

  !> Diagonalize a real tridiagonal symmetric matrix \( \mathbf{T} \) such that
  !> \[
  !>    \mathbf{T} \mathbf{x_i} = \lambda_i \mathbf{x_i},
  !> \]
  !> where \( \lambda_i \) is the \(i\)'th eigenvalue and \( \mathbf{x_i} \) \(i\)'th eigenvector.
  !>
  !> This routine acts on the arrays as expected by the LAPACK routine [[dstedc]].
  subroutine dstedc_wrapper(compz, N, D, E, Z, lwork_in, liwork_in)
    !> Determine what is to be computed
    !> 'N':  Compute eigenvalues only.
    !> 'I':  Compute eigenvectors of tridiagonal matrix also.
    !> 'V':  Compute eigenvectors of original dense symmetric
    !>       matrix also.  On entry, eigenvectors contains the orthogonal
    !>       matrix used to reduce the original matrix to
    !>       tridiagonal form.
    character(len=1), intent(in) :: compz
    !> Length of the diagonal
    integer, intent(in) :: N
    !> On entry, the \( N \) diagonal elements of the tridiagonal matrix \( \mathbf{T} \).
    !> On exit, the eigenvalues in ascending order
    real(dp), intent(inout) :: D(N)
    !> On entry, the \( N - 1 \) subdiagonal elements of the tridiagonal matrix \( \mathbf{T} \).
    !> On exit, it has been destroyed.
    real(dp), intent(inout) :: E(N - 1)
    !> The content `Z` has on entry and exit depends on `compz`.
    !>
    !> - If `compz = N`, `Z` is not referenced
    !>
    !> - If `compz = I`, on exit, `Z` contains the orthonormal eigenvectors of
    !>   the symmetric tridiagonal matrix \( \mathbf{T} \).
    !>
    !> - If `compz = V`,
    !>
    !>     - on entry, `Z` contains the orthogonal matrix used in the reduction to the tridiagonal form,
    !>
    !>     - on exit, `Z` contains the orthonormal eigenvectors of the orginal matrix.
    real(dp), intent(inout) :: Z(N, N)
    !> Dimension of the work array WORK. If not present, a work space query is done.
    integer, intent(in), optional :: lwork_in
    !> Dimension of the work array IWORK. If not present, a work space query is done.
    integer, intent(in), optional :: liwork_in

    ! work space for dstedc
    integer :: lwork, liwork, info
    real(dp), allocatable :: work(:)
    integer, allocatable :: iwork(:)

    call assert(any(compz == allowed_compz), 'compz is not one of the allowed chracters ("N", "I", "V").')

    if (present(lwork_in) .and. present(liwork_in)) then
      lwork = lwork_in
      liwork = liwork_in
    else
      call xstedc_wsq(compz, N, D, E, Z, lwork, liwork)
    end if

    allocate(work(lwork))
    allocate(iwork(liwork))
    call dstedc(compz, N, D, E, Z, N, work, lwork, iwork, liwork, info)
    call terminate_if_false(info == 0, 'dstedc failed.')
  end subroutine dstedc_wrapper

  !> Diagonalize a real tridiagonal symmetric matrix \( \mathbf{T} \) such that
  !> \[
  !>    \mathbf{T} \mathbf{x_i} = \lambda_i \mathbf{x_i},
  !> \]
  !> where \( \lambda_i \) is the \(i\)'th eigenvalue and \( \mathbf{x_i} \) \(i\)'th eigenvector.
  !>
  !> This routine acts on the arrays as expected by the LAPACK routine [[dstedc]].
  subroutine dstedc_wrapper_workspace_query(compz, N, D, E, Z, lwork, liwork)
    !> Determine what is to be computed
    !> 'N':  Compute eigenvalues only.
    !> 'I':  Compute eigenvectors of tridiagonal matrix also.
    !> 'V':  Compute eigenvectors of original dense symmetric
    !>       matrix also.  On entry, `eigenvectors` contains the orthogonal
    !>       matrix used to reduce the original matrix to
    !>       tridiagonal form.
    character(len=1), intent(in) :: compz
    !> Length of the diagonal
    integer, intent(in) :: N
    !> On entry, the \( N \) diagonal elements of the tridiagonal matrix \( \mathbf{T} \).
    !> On exit, the eigenvalues in ascending order
    real(dp), intent(inout) :: D(N)
    !> On entry, the \( N - 1 \) subdiagonal elements of the tridiagonal matrix \( \mathbf{T} \).
    !> On exit, `E` has been destroyed.
    real(dp), intent(inout) :: E(N - 1)
    !> The content `Z` has on entry and exit depends on `compz`.
    !>
    !> - If `compz = N`, `Z` is not referenced
    !>
    !> - If `compz = I`, on exit, `Z` contains the orthonormal eigenvectors of
    !>   the symmetric tridiagonal matrix \( \mathbf{T} \).
    !>
    !> - If `compz = V`,
    !>
    !>     - on entry, `Z` contains the orthogonal matrix used in the reduction to the tridiagonal form,
    !>
    !>     - on exit, `Z` contains the orthonormal eigenvectors of the orginal matrix.
    real(dp), intent(inout) :: Z(N, N)
    !> Dimension of the work array WORK
    integer, intent(out) :: lwork
    !> Dimension of the work array IWORK
    integer, intent(out) :: liwork

    ! work space for dstedc
    integer :: info, iwork(1)
    real(dp) :: work(1)

    call assert(any(compz == allowed_compz), 'compz is not one of the allowed chracters ("N", "I", "V").')

    lwork = -1
    liwork = -1
    call dstedc(compz, N, D, E, Z, N, work, lwork, iwork, liwork, info)
    call terminate_if_false(info == 0, 'dstedc failed for workspace query.')

    lwork = nint(work(1))
    liwork = iwork(1)
  end subroutine dstedc_wrapper_workspace_query

end module diagonalize_tridiagonal
