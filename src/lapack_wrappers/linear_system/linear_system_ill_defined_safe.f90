!> Module with LAPACK wrappers for solving system of linear equations
!>    \[ A * x = y \]
!> where `A` is a hermitian positive definite matrix
module linear_system_ill_defined_safe
  use asserts, only: assert
  use modmpi, only: terminate_if_false
  use lapack_f95_interfaces, only: zgelsd
  use lapack_workspaces, only: lapack_workspace_complex_dp_t
  use precision, only: dp, i32

  implicit none

  private

  public :: ill_defined_safe_solve

  interface ill_defined_safe_solve
    module procedure :: ill_defined_safe_solve_complex_dp
  end interface ill_defined_safe_solve

  real(dp), parameter :: threshold_default = 1.0e-6_dp

contains

  !> Solve the system of linear equations \(A \cdot x = y\) for complex double precision variables,
  !> where \(A\) is a general (not necessarily square or full-rank) complex matrix.
  !>
  !> This subroutine uses a singular value decomposition (SVD)-based method to compute
  !> the Moore–Penrose pseudoinverse of A, solving the system in a least-squares sense.
  !> Singular values below a certain threshold are discarded during the inversion.
  !>
  subroutine ill_defined_safe_solve_complex_dp(A, y, workspace, threshold)
    !> Matrix \(A\) of the \(A \cdot x = y\) problem
    complex(dp), contiguous, intent(in) :: A(:, :)
    !> On entry, the (\y\) matrix, on exit the x of the \(A \cdot x = y\) problem
    complex(dp), contiguous, intent(inout) :: y(:, :)
    !> workspace (optional)
    type(lapack_workspace_complex_dp_t), target, optional, intent(inout) :: workspace
    !> The threshold to consider a value 0 with respect to the maximum singular value
    real(dp), intent(in), optional :: threshold

    ! Local version of the optionals
    type(lapack_workspace_complex_dp_t), pointer :: workspace_fptr
    type(lapack_workspace_complex_dp_t), target  :: workspace_local
    real(dp)                                     :: threshold_local


    integer(i32) :: m, n, nhrs, rank, info
    real(dp), allocatable :: singular_values(:)
    complex(dp), allocatable :: A_copy(:, :)

    m    = size(A,1)
    n    = size(A,2)
    nhrs = size(y,2)
    A_copy = A

    allocate(singular_values(min(m,n)))

    ! Associate the workspace to either a local or an external workspace
    if (present(workspace)) then
      workspace_fptr => workspace
    else
      workspace_fptr => workspace_local
    end if

    if (present(threshold)) then
      threshold_local = threshold
    else
      threshold_local = threshold_default
    end if

    ! Get optimal workspace if no workspace is present or has not been inited
    if (.not. workspace_local%computed()) then
      call workspace_local%initialize(lrwork=1, liwork=1, lwork=1)
      call workspace_local%allocate_workspace()
      call zgelsd(m, n, nhrs, A_copy, m, y, max(m,n), singular_values, threshold_local, rank, workspace_local%work, -1, &
                  workspace_local%rwork, workspace_local%iwork, info)
      call terminate_if_false(info == 0, "ill_defined_safe_solve_complex_dp: subspace query failed")
      call workspace_local%reset(.true.)
    end if

    ! Solve
    call zgelsd(m, n, nhrs, A_copy, m, y, max(m,n), singular_values, threshold_local, rank, workspace_local%work, size(workspace_local%work), &
                workspace_local%rwork, workspace_local%iwork, info)

    call terminate_if_false(info == 0, "ill_defined_safe_solve_complex_dp: zgelsd failed")

    nullify(workspace_fptr)

  end subroutine ill_defined_safe_solve_complex_dp
end module linear_system_ill_defined_safe
